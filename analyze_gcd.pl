#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $dir = '.';
my @files;

my $mode = 'allele';
my $help = 0;
my $debug = 0;
GetOptions(
    'dir=s' => \$dir,
    'mode=s' => \$mode,
    'debug!' => \$debug,
    'help|h' => \$help,
) or die "Usage: $0 [--dir <directory>] [--mode allele|serology] [--debug] [--help]\n";

@files = glob("$dir/ION-*-[DC].csv");

if ($help) {
    print "Usage: $0 [--dir <directory>] [--mode allele|serology] [--debug] [--help]\n";
    print "\n";
    print "Analyzes HLA typing data from ION CSV files to determine typing resolution quality.\n";
    print "\n";
    print "Options:\n";
    print "  --dir <dir>     Input directory with ION-*-[DC].csv files (default: current directory)\n";
    print "  --mode allele   Require allele-level resolution (default)\n";
    print "  --mode serology Count donors with serology-only loci (no DNA)\n";
    print "  --debug         Show debug info for each donor\n";
    print "  --help          Show this help\n";
    print "\n";
    print "Example:\n";
    print "  $0 --dir examples/input --mode allele\n";
    exit 0;
}

$mode = lc $mode;
die "Invalid --mode '$mode'. Use 'allele' or 'serology'.\n" unless $mode eq 'allele' || $mode eq 'serology';

my $output_dir = "$dir/analysis_results";

my $missing_label = ($mode eq 'allele') ? 'No Allele (SER)' : 'No DNA (SER)';
my $csv_missing_label = ($mode eq 'allele') ? 'all_allele_high' : 'No_DNA_SER';

print "Mode: $mode\n";

die "No CSV files found!" unless @files;

# Create output directory
mkdir $output_dir unless -d $output_dir;

# Fields that must have allele-level resolution
my @required_fields = qw(DNA_A1 DNA_A2 DNA_B1 DNA_B2 DNA_C1 DNA_C2 DRB11 DRB12 DQB11 DQB12);

# IONs to exclude from analysis
# 8405 (KR) is suspected to reporting serology-only data for all donors derived from DNA-typed samples
my %exclude_ions = map { $_ => 1 } qw(8405);

my %total_stats;  # Overall statistics

# Current date for age calculation
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $today_year = $year + 1900;
my $today_month = $mon + 1;
my $today_mday = $mday;

foreach my $file (sort @files) {
    print "Processing: $file\n";
    
    my ($ion) = $file =~ /ION-(\d+)-([DC])\.csv/;
    my $type = $2;
    
    # Skip excluded IONs
    if ($exclude_ions{$ion}) {
        print "  Skipping excluded ION: $ion\n";
        next;
    }
    
    open(my $fh, '<', $file) or do {
        warn "Cannot open $file: $!";
        next;
    };
    
    # Read header
    my $header = <$fh>;
    chomp $header;
    
    # Split header more robustly
    my @columns = split(',', $header);
    
    # Find column indices only for fields we need
    my %col_idx;
    my @serologic_fields = qw(A1 A2 B1 B2 C1 C2 DR1 DR2 DQ1 DQ2);
    my @all_lookup_fields = ('BIRTH_DATE', 'ID', @required_fields, @serologic_fields);
    for my $i (0 .. $#columns) {
        $col_idx{$columns[$i]} = $i if grep { $_ eq $columns[$i] } @all_lookup_fields;
    }
    
    # Verify that all required fields exist
    my $fields_missing = 0;
    foreach my $field (@required_fields) {
        unless (exists $col_idx{$field}) {
            warn "Field '$field' not found in $file (skipping this file)";
            $fields_missing = 1;
        }
    }
    
    my $birth_date_idx = $col_idx{BIRTH_DATE} // -1;
    if ($birth_date_idx < 0) {
        warn "BIRTH_DATE not found in $file (skipping this file)";
        $fields_missing = 1;
    }
    
    if ($fields_missing) {
        close($fh);
        next;
    }
    
    # Statistics per ION/Type/Age
    my %local_stats;
    my $line_count = 0;
    
    while (my $line = <$fh>) {
        chomp $line;
        next if !$line;
        
        $line_count++;
        
        # Fast field extraction - only split, don't create full array
        my @fields = split(',', $line, -1);  # -1 keeps trailing empty fields
        
        # Skip if not enough fields
        next if @fields < (scalar(keys %col_idx) + 10);
        
        # Calculate birth date and age
        my $birth_date_str = $fields[$birth_date_idx] // '';
        my $age = '';
        
        # Debug: Track why dates fail to parse
        if (!$birth_date_str) {
            # Empty field
        } elsif ($birth_date_str =~ /^(\d{4})-(\d{2})-(\d{2})$/) {
            my $birth_year = $1;
            my $birth_month = $2;
            my $birth_mday = $3;
            
            $age = $today_year - $birth_year;
            # Correct for birthday not yet passed
            if ($birth_month > $today_month || 
                ($birth_month == $today_month && $birth_mday > $today_mday)) {
                $age--;
            }
        } else {
            # Date format doesn't match YYYY-MM-DD - log example
            # (we'll count these separately)
        }
        
        $age = 'unknown' if $age eq '';
        
        my $has_full_resolution = 1;

        if ($mode eq 'serology') {
            $has_full_resolution = 1;
            my %loci_ser = (
                A  => ['A1', 'A2'],
                B  => ['B1', 'B2'],
                C  => ['C1', 'C2'],
                DR => ['DR1', 'DR2'],
                DQ => ['DQ1', 'DQ2'],
            );
            my %loci_dna = (
                A  => ['DNA_A1', 'DNA_A2'],
                B  => ['DNA_B1', 'DNA_B2'],
                C  => ['DNA_C1', 'DNA_C2'],
                DR => ['DRB11',  'DRB12'],
                DQ => ['DQB11',  'DQB12'],
            );
            foreach my $locus (keys %loci_ser) {
                my ($ser_field_1, $ser_field_2) = @{$loci_ser{$locus}};
                my ($dna_field_1, $dna_field_2) = @{$loci_dna{$locus}};
                my $ser_val_1 = $fields[$col_idx{$ser_field_1}] // '';
                my $dna_val_1 = $fields[$col_idx{$dna_field_1}] // '';
                my $dna_val_2 = $fields[$col_idx{$dna_field_2}] // '';
                # Has serologic assignment at this locus but no DNA assignment -> serology-only
                if ($ser_val_1 ne '' && $dna_val_1 eq '' && $dna_val_2 eq '') {
                    $has_full_resolution = 0;
                    last;
                }
            }
        } else {
            my %loci = (
                DNA_A => ['DNA_A1', 'DNA_A2'],
                DNA_B => ['DNA_B1', 'DNA_B2'],
                DNA_C => ['DNA_C1', 'DNA_C2'],
                DRB1  => ['DRB11',  'DRB12'],
                DQB1  => ['DQB11',  'DQB12'],
            );

            foreach my $locus (keys %loci) {
                my ($field_1, $field_2) = @{$loci{$locus}};
                my $val_1 = $fields[$col_idx{$field_1}] // '';
                my $val_2 = $fields[$col_idx{$field_2}] // '';

                # First allele must be present and in allele-level resolution
                if (!is_allele_level($val_1)) {
                    $has_full_resolution = 0;
                    last;
                }

                # Second allele can be empty only if first is allele-level (homozygous)
                if ($val_2 ne '' && !is_allele_level($val_2)) {
                    $has_full_resolution = 0;
                    last;
                }
            }
            
            # Debug output for allele mode
            if ($debug && !$has_full_resolution) {
                my $donor_id = $fields[$col_idx{ID}] // '(no ID)';
                debug_allele_donor($donor_id, $age, \@fields, \%col_idx);
            }
        }

        # Count
        $local_stats{$age}{'all_rows'}++;
        $local_stats{$age}{'missing'}++ unless $has_full_resolution;

        $total_stats{$type}{$age}{'all_rows'}++;
        $total_stats{$type}{$age}{'missing'}++ unless $has_full_resolution;
        
        # Debug output for serology mode
        if ($debug && $mode eq 'serology' && !$has_full_resolution) {
            my $donor_id = $fields[$col_idx{ID}] // '(no ID)';
            debug_serology_donor($donor_id, $age, \@fields, \%col_idx);
        }
    }
    
    close($fh);
    
    # Output results for this file
    printf "%-10s %15s %15s %12s\n", "Age", "Total", $missing_label, "Percent";
    print "-" x 80 . "\n";
    
    my $total_all = 0;
    my $total_missing = 0;
    
    foreach my $age (sort { 
        my $a_num = ($a eq 'unknown') ? 999 : $a;
        my $b_num = ($b eq 'unknown') ? 999 : $b;
        $a_num <=> $b_num
    } keys %local_stats) {
        my $all = $local_stats{$age}{'all_rows'} // 0;
        my $missing = $local_stats{$age}{'missing'} // 0;
        my $percent = $all > 0 ? sprintf("%.1f", $missing * 100 / $all) : 0;
        printf "%-10s %15d %15d %11s%%\n", $age, $all, $missing, $percent;
        $total_all += $all;
        $total_missing += $missing;
    }
    
    print "-" x 80 . "\n";
    my $percent = $total_all > 0 ? sprintf("%.1f", $total_missing * 100 / $total_all) : 0;
    printf "%-10s %15d %15d %11s%%\n", "TOTAL", $total_all, $total_missing, $percent;
    
    # Write CSV file for this ION
    write_csv_file("$output_dir/ION-$ion-$type.csv", \%local_stats, $csv_missing_label);
}

# Summary statistics
print "\n\n";
print "=" x 80 . "\n";
print "SUMMARY (D-Donors)\n";
print "=" x 80 . "\n";
printf "%-10s %15s %15s %12s\n", "Age", "Total", $missing_label, "Percent";
print "-" x 80 . "\n";

my (@summary_d, @summary_c);
my $grand_total_all = 0;
my $grand_total_missing = 0;

if (exists $total_stats{'D'}) {
    foreach my $age (sort { 
        my $a_num = ($a eq 'unknown') ? 999 : $a;
        my $b_num = ($b eq 'unknown') ? 999 : $b;
        $a_num <=> $b_num
    } keys %{$total_stats{'D'}}) {
        my $all = $total_stats{'D'}{$age}{'all_rows'} // 0;
        my $missing = $total_stats{'D'}{$age}{'missing'} // 0;
        my $percent = $all > 0 ? sprintf("%.1f", $missing * 100 / $all) : 0;
        printf "%-10s %15d %15d %11s%%\n", $age, $all, $missing, $percent;
        push @summary_d, [$age, $all, $missing, $percent];
        $grand_total_all += $all;
        $grand_total_missing += $missing;
    }
}

print "-" x 80 . "\n";
my $percent = $grand_total_all > 0 ? sprintf("%.1f", $grand_total_missing * 100 / $grand_total_all) : 0;
printf "%-10s %15d %15d %11s%%\n", "TOTAL", $grand_total_all, $grand_total_missing, $percent;
push @summary_d, ["TOTAL", $grand_total_all, $grand_total_missing, $percent];

# Write CSV file for D summary
write_summary_csv("$output_dir/summary_D-Donors.csv", \@summary_d, $csv_missing_label);

print "\n";
print "=" x 80 . "\n";
print "SUMMARY (C-Donors)\n";
print "=" x 80 . "\n";
printf "%-10s %15s %15s %12s\n", "Age", "Total", $missing_label, "Percent";
print "-" x 80 . "\n";

$grand_total_all = 0;
$grand_total_missing = 0;

if (exists $total_stats{'C'}) {
    foreach my $age (sort { 
        my $a_num = ($a eq 'unknown') ? 999 : $a;
        my $b_num = ($b eq 'unknown') ? 999 : $b;
        $a_num <=> $b_num
    } keys %{$total_stats{'C'}}) {
        my $all = $total_stats{'C'}{$age}{'all_rows'} // 0;
        my $missing = $total_stats{'C'}{$age}{'missing'} // 0;
        my $percent = $all > 0 ? sprintf("%.1f", $missing * 100 / $all) : 0;
        printf "%-10s %15d %15d %11s%%\n", $age, $all, $missing, $percent;
        push @summary_c, [$age, $all, $missing, $percent];
        $grand_total_all += $all;
        $grand_total_missing += $missing;
    }
}

print "-" x 80 . "\n";
$percent = $grand_total_all > 0 ? sprintf("%.1f", $grand_total_missing * 100 / $grand_total_all) : 0;
printf "%-10s %15d %15d %11s%%\n", "TOTAL", $grand_total_all, $grand_total_missing, $percent;
push @summary_c, ["TOTAL", $grand_total_all, $grand_total_missing, $percent];

# Write CSV file for C summary
write_summary_csv("$output_dir/summary_C-Donors.csv", \@summary_c, $csv_missing_label);

# Overall summary
print "\n";
print "=" x 80 . "\n";
print "SUMMARY (D + C Combined)\n";
print "=" x 80 . "\n";
printf "%-10s %15s %15s %12s\n", "Age", "Total", $missing_label, "Percent";
print "-" x 80 . "\n";

my %combined_stats;
foreach my $type (keys %total_stats) {
    foreach my $age (keys %{$total_stats{$type}}) {
        $combined_stats{$age}{'all_rows'} += $total_stats{$type}{$age}{'all_rows'} // 0;
        $combined_stats{$age}{'missing'} += $total_stats{$type}{$age}{'missing'} // 0;
    }
}

$grand_total_all = 0;
$grand_total_missing = 0;
my @summary_combined;

foreach my $age (sort { 
    my $a_num = ($a eq 'unknown') ? 999 : $a;
    my $b_num = ($b eq 'unknown') ? 999 : $b;
    $a_num <=> $b_num
} keys %combined_stats) {
    my $all = $combined_stats{$age}{'all_rows'} // 0;
    my $missing = $combined_stats{$age}{'missing'} // 0;
    my $percent = $all > 0 ? sprintf("%.1f", $missing * 100 / $all) : 0;
    printf "%-10s %15d %15d %11s%%\n", $age, $all, $missing, $percent;
    push @summary_combined, [$age, $all, $missing, $percent];
    $grand_total_all += $all;
    $grand_total_missing += $missing;
}

print "-" x 80 . "\n";
$percent = $grand_total_all > 0 ? sprintf("%.1f", $grand_total_missing * 100 / $grand_total_all) : 0;
printf "%-10s %15d %15d %11s%%\n", "TOTAL", $grand_total_all, $grand_total_missing, $percent;
push @summary_combined, ["TOTAL", $grand_total_all, $grand_total_missing, $percent];

# Write CSV file for combined summary
write_summary_csv("$output_dir/summary_Overall.csv", \@summary_combined, $csv_missing_label);

print "\nAnalysis completed.\n";
print "CSV files were saved in $output_dir\n";

# ============================================================================
# Helper functions
# ============================================================================

sub write_csv_file {
    my ($filename, $stats_ref, $missing_label) = @_;
    
    open(my $csv_fh, '>', $filename) or do {
        warn "Cannot write CSV file $filename: $!";
        return;
    };
    
    # Write header
    print $csv_fh "Age,Total,$missing_label,Percent\n";
    
    # Write data
    foreach my $age (sort { 
        my $a_num = ($a eq 'unknown') ? 999 : $a;
        my $b_num = ($b eq 'unknown') ? 999 : $b;
        $a_num <=> $b_num
    } keys %{$stats_ref}) {
        my $all = $stats_ref->{$age}{'all_rows'} // 0;
        my $missing = $stats_ref->{$age}{'missing'} // 0;
        my $percent = $all > 0 ? sprintf("%.1f", $missing * 100 / $all) : 0;
        print $csv_fh "$age,$all,$missing,$percent\n";
    }
    
    close($csv_fh);
    print "  Written: $filename\n";
}

sub write_summary_csv {
    my ($filename, $data_array, $missing_label) = @_;
    
    open(my $csv_fh, '>', $filename) or do {
        warn "Cannot write CSV file $filename: $!";
        return;
    };
    
    # Write header
    print $csv_fh "Age,Total,$missing_label,Percent\n";
    
    # Write data
    foreach my $row (@{$data_array}) {
        my ($age, $all, $missing, $percent) = @{$row};
        print $csv_fh "$age,$all,$missing,$percent\n";
    }
    
    close($csv_fh);
    print "  Written: $filename\n";
}

sub is_allele_level {
    my ($value) = @_;

    return 0 unless defined $value;
    return 0 if $value eq '';

    # 2-4 colon-separated numeric fields, optional trailing letter at the end
    return ($value =~ /^\d+(?::\d+){1,3}[A-Z]?$/) ? 1 : 0;
}

sub debug_serology_donor {
    my ($donor_id, $age, $fields_ref, $col_idx_ref) = @_;
    my @fields = @$fields_ref;
    my %col_idx = %$col_idx_ref;
    
    my %loci_ser = (
        A  => ['A1', 'A2'],
        B  => ['B1', 'B2'],
        C  => ['C1', 'C2'],
        DR => ['DR1', 'DR2'],
        DQ => ['DQ1', 'DQ2'],
    );
    my %loci_dna = (
        A  => ['DNA_A1', 'DNA_A2'],
        B  => ['DNA_B1', 'DNA_B2'],
        C  => ['DNA_C1', 'DNA_C2'],
        DR => ['DRB11',  'DRB12'],
        DQ => ['DQB11',  'DQB12'],
    );
    
    print "  $donor_id (age $age): serology-only loci\n";
    foreach my $locus (sort keys %loci_ser) {
        my ($ser_field_1, $ser_field_2) = @{$loci_ser{$locus}};
        my ($dna_field_1, $dna_field_2) = @{$loci_dna{$locus}};
        my $ser_val_1 = $fields[$col_idx{$ser_field_1}] // '';
        my $ser_val_2 = $fields[$col_idx{$ser_field_2}] // '';
        my $dna_val_1 = $fields[$col_idx{$dna_field_1}] // '';
        my $dna_val_2 = $fields[$col_idx{$dna_field_2}] // '';
        
        my $ser_display = ($ser_val_1 ne '' || $ser_val_2 ne '') ? 
            ($ser_val_1 ne '' ? $ser_val_1 : '-') . ',' . ($ser_val_2 ne '' ? $ser_val_2 : '-') :
            '-,-';
        my $dna_display = ($dna_val_1 ne '' || $dna_val_2 ne '') ? 
            ($dna_val_1 ne '' ? $dna_val_1 : '-') . ',' . ($dna_val_2 ne '' ? $dna_val_2 : '-') :
            '-,-';
        
        my $marker = ($ser_val_1 ne '' && $dna_val_1 eq '' && $dna_val_2 eq '') ? ' [SERO-ONLY]' : '';
        printf "    %-3s Ser=%-15s DNA=%s%s\n", $locus . ':', $ser_display, $dna_display, $marker;
    }
}

sub debug_allele_donor {
    my ($donor_id, $age, $fields_ref, $col_idx_ref) = @_;
    my @fields = @$fields_ref;
    my %col_idx = %$col_idx_ref;
    
    my %loci = (
        DNA_A => ['DNA_A1', 'DNA_A2'],
        DNA_B => ['DNA_B1', 'DNA_B2'],
        DNA_C => ['DNA_C1', 'DNA_C2'],
        DRB1  => ['DRB11',  'DRB12'],
        DQB1  => ['DQB11',  'DQB12'],
    );
    
    print "  $donor_id (age $age): allele-level resolution lacking\n";
    foreach my $locus (sort keys %loci) {
        my ($field_1, $field_2) = @{$loci{$locus}};
        my $val_1 = $fields[$col_idx{$field_1}] // '';
        my $val_2 = $fields[$col_idx{$field_2}] // '';
        my $val_1_ok = is_allele_level($val_1);
        my $val_2_ok = ($val_2 eq '') || is_allele_level($val_2);
        
        my $display_1 = ($val_1 ne '') ? $val_1 : '-';
        my $display_2 = ($val_2 ne '') ? $val_2 : '-';
        my $reason = '';
        if (!$val_1_ok) {
            $reason = " [1st invalid: '$display_1']";
        } elsif (!$val_2_ok) {
            $reason = " [2nd invalid: '$display_2']";
        }
        printf "    %-5s %s, %s%s\n", $locus . ':', $display_1, $display_2, $reason;
    }
}
