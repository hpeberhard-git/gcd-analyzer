#!/usr/bin/perl -w

use strict;
use warnings;
use Text::CSV_XS;
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
    print "  --mode allele   Count donors with allele-level resolution (default)\n";
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

our %loci_ser = (
    A  => ['A1', 'A2'],
    B  => ['B1', 'B2'],
    C  => ['C1', 'C2'],
    DR => ['DR1', 'DR2'],
    DQ => ['DQ1', 'DQ2'],
);

our %loci_dna = (
    A  => ['DNA_A1', 'DNA_A2'],
    B  => ['DNA_B1', 'DNA_B2'],
    C  => ['DNA_C1', 'DNA_C2'],
    DR => ['DRB11',  'DRB12'],
    DQ => ['DQB11',  'DQB12'],
);

# Fields that must have allele-level resolution
our @locus_list = qw(A B C DR DQ);

my @ser_fields = map { @{ $loci_ser{$_} } } @locus_list;
my @dna_fields = map { @{ $loci_dna{$_} } } @locus_list;
my @hla_fields = (@dna_fields, @ser_fields);

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

    my $csv = Text::CSV_XS->new({
       binary      => 1,
       auto_diag   => 0,
       sep_char    => ',',
       quote_char  => '"',
       escape_char => '"',
    });

    open(my $fh, '<', $file) or do {
        warn "Cannot open $file: $!";
        next;
    };

    # Read header
    my $header = <$fh>;
    $header =~ s/\s+$//;

    if ($header =~ /^DELETE ALL/) {
        warn "Skip 'empty' $file\n";
        close($fh);
        next;
    }

    # split header robustly
    $csv->parse($header),
    my @columns = $csv->fields();
    my $col_count = @columns;

    # Find column indices only for fields we need
    my %col_idx;
    my @all_lookup_fields = ('BIRTH_DATE', 'ID', @hla_fields);
    for my $i (0 .. $#columns) {
        $col_idx{$columns[$i]} = $i if grep { $_ eq $columns[$i] } @all_lookup_fields;
    }

    # Verify that all hla fields exist
    my $fields_missing = 0;
    foreach my $field (@hla_fields) {
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
        $line =~ s/\s+$//;
        next if !$line;

        $line_count++;

        # split line robustly
        $csv->parse($line);
        my @fields = $csv->fields();
        my $field_count = @fields;

        # Skip line if structure is unexpected
        if ($field_count != $col_count) {
           warn "line $line_count skipped: expected $col_count fields but found $field_count\n";
           next;
        }

        # Calculate donor's age
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
            warn "line $line_count: unexpected format of date of birth: $birth_date_str\n";
        }

        $age = 'unknown' if $age eq '';

        my $has_full_resolution = 1;

        if ($mode eq 'serology') {
            $has_full_resolution = 1;

            foreach my $locus (@locus_list) {
                my ($ser_field_1, $ser_field_2) = @{$loci_ser{$locus}};
                my ($dna_field_1, $dna_field_2) = @{$loci_dna{$locus}};

                my $ser_val_1 = &normalize_field($fields[ $col_idx{$ser_field_1} ]);
                my $ser_val_2 = &normalize_field($fields[ $col_idx{$ser_field_2} ]);
                my $dna_val_1 = &normalize_field($fields[ $col_idx{$dna_field_1} ]);
                my $dna_val_2 = &normalize_field($fields[ $col_idx{$dna_field_2} ]);

                my $has_ser = ($ser_val_1 ne '' || $ser_val_2 ne '');
                my $has_dna = ($dna_val_1 ne '' || $dna_val_2 ne '');

                # Has serology at this locus but no DNA assignment -> not full resolution
                if ($has_ser && !$has_dna) {
                   $has_full_resolution = 0;
                   last;
                }
            }
        }
        else {
            foreach my $locus (@locus_list) {
                my ($dna_field_1, $dna_field_2) = @{$loci_dna{$locus}};

                my $dna_val_1 = &normalize_field($fields[$col_idx{$dna_field_1}]);
                my $dna_val_2 = &normalize_field($fields[$col_idx{$dna_field_2}]);

                # First allele must be present and in allele-level resolution
                if (!&is_allele_level($dna_val_1)) {
                    $has_full_resolution = 0;
                    last;
                }

                # Second allele can be empty only if first is allele-level (homozygous)
                if ($dna_val_2 ne '' && !is_allele_level($dna_val_2)) {
                    $has_full_resolution = 0;
                    last;
                }
            }
        }

        # Debug output
        if ($debug && !$has_full_resolution) {
           my $donor_id = $fields[$col_idx{'ID'}] // '(no ID)';
           ($mode eq 'serology')
              ? &debug_serology_donor($donor_id, $age, \@fields, \%col_idx)
              : &debug_allele_donor($donor_id, $age, \@fields, \%col_idx);
        }

        # Count
        $local_stats{$age}{'all_rows'}++;
        $local_stats{$age}{'missing'}++ unless $has_full_resolution;

        $total_stats{$type}{$age}{'all_rows'}++;
        $total_stats{$type}{$age}{'missing'}++ unless $has_full_resolution;
    }

    close($fh);

    # Output results for this file
    my $rows_ref = build_stats_rows(\%local_stats);

    &print_stats_table(
        rows_ref      => $rows_ref,
        missing_label => $missing_label,
    );

    # Write CSV file for this ION (without TOTAL row)
    &write_csv_rows("$output_dir/ION-$ion-$type.csv", $rows_ref, $csv_missing_label, 0);
}

# Summary statistics
print "\n\n";

# D summary
&output_summary_report(
    title            => "SUMMARY (D-Donors)",
    stats_ref        => (exists $total_stats{'D'} ? $total_stats{'D'} : {}),
    missing_label    => $missing_label,
    csv_file         => "$output_dir/summary_D-Donors.csv",
    csv_missing_label=> $csv_missing_label,
);

print "\n";

# C summary
&output_summary_report(
    title            => "SUMMARY (C-Donors)",
    stats_ref        => (exists $total_stats{'C'} ? $total_stats{'C'} : {}),
    missing_label    => $missing_label,
    csv_file         => "$output_dir/summary_C-Donors.csv",
    csv_missing_label=> $csv_missing_label,
);

print "\n";

# Overall summary (D + C combined)
my %combined_stats;
foreach my $type (keys %total_stats) {
    foreach my $age (keys %{$total_stats{$type}}) {
        $combined_stats{$age}{'all_rows'} += $total_stats{$type}{$age}{'all_rows'} // 0;
        $combined_stats{$age}{'missing'}  += $total_stats{$type}{$age}{'missing'}  // 0;
    }
}

&output_summary_report(
    title            => "SUMMARY (D + C Combined)",
    stats_ref        => \%combined_stats,
    missing_label    => $missing_label,
    csv_file         => "$output_dir/summary_Overall.csv",
    csv_missing_label=> $csv_missing_label,
);

print "\nAnalysis completed.\n";
print "CSV files were saved in $output_dir\n";

# ============================================================================
# Helper functions
# ============================================================================


sub normalize_field {
    my ($value) = @_;

    $value //= '';

    # Trim outer whitespace
    $value =~ s/^\s+|\s+$//g;

    # Remove one pair of surrounding quotes
    $value =~ s/^"(.*)"$/$1/;
    $value =~ s/^'(.*)'$/$1/;

    # Trim again after unquoting
    $value =~ s/^\s+|\s+$//g;

    return $value;
}

sub sort_age_keys {
    my ($stats_ref) = @_;
    return sort {
        ($a eq 'unknown') <=> ($b eq 'unknown')
            ||
        $a <=> $b
    } keys %{$stats_ref};
}

sub build_stats_rows {
    my ($stats_ref) = @_;

    my $total_all     = 0;
    my $total_missing = 0;
    my @rows;

    foreach my $age (sort_age_keys($stats_ref)) {
        my $all     = $stats_ref->{$age}{'all_rows'} // 0;
        my $missing = $stats_ref->{$age}{'missing'}  // 0;
        my $percent = $all > 0 ? sprintf("%.1f", $missing * 100 / $all) : 0;

        push @rows, [$age, $all, $missing, $percent];
        $total_all     += $all;
        $total_missing += $missing;
    }

    my $total_percent = $total_all > 0 ? sprintf("%.1f", $total_missing * 100 / $total_all) : 0;
    push @rows, ["TOTAL", $total_all, $total_missing, $total_percent];

    return \@rows;
}

sub print_stats_table {
    my (%args) = @_;
    my $rows_ref      = $args{rows_ref}      // [];
    my $missing_label = $args{missing_label} // 'Missing';
    my $title         = $args{title};                 # optional
    my $banner        = $args{banner} // 0;           # 1 => print === lines

    if (defined $title) {
        if ($banner) {
            print "=" x 80 . "\n";
            print "$title\n";
            print "=" x 80 . "\n";
        } else {
            print "$title\n";
        }
    }

    printf "%-10s %15s %15s %12s\n", "Age", "Total", $missing_label, "Percent";
    print "-" x 80 . "\n";

    for my $row (@{$rows_ref}) {
        my ($age, $all, $missing, $percent) = @{$row};

        if ($age eq 'TOTAL') {
            print "-" x 80 . "\n";
        }

        printf "%-10s %15d %15d %11s%%\n", $age, $all, $missing, $percent;
    }
}

sub write_csv_rows {
    my ($filename, $rows_ref, $missing_label, $include_total) = @_;

    # Default: include TOTAL row
    $include_total = 1 if !defined $include_total;

    open(my $csv_fh, '>', $filename) or do {
        warn "Cannot write CSV file $filename: $!";
        return;
    };

    # Write header
    print $csv_fh "Age,Total,$missing_label,Percent\n";

    # Write data
    foreach my $row (@{$rows_ref}) {
        my ($age, $all, $missing, $percent) = @{$row};

        next if (!$include_total && $age eq 'TOTAL');

        print $csv_fh "$age,$all,$missing,$percent\n";
    }

    close($csv_fh);
    print "  Written: $filename\n";
}

sub print_summary_section {
    my (%args) = @_;
    my $title         = $args{title} // 'SUMMARY';
    my $stats_ref     = $args{stats_ref} // {};
    my $missing_label = $args{missing_label} // 'Missing';

    my $rows_ref = build_stats_rows($stats_ref);

    print_stats_table(
        title         => $title,
        banner        => 1,
        rows_ref      => $rows_ref,
        missing_label => $missing_label,
    );

    return $rows_ref;
}

sub output_summary_report {
    my (%args) = @_;
    my $title            = $args{title} // 'SUMMARY';
    my $stats_ref        = $args{stats_ref} // {};
    my $missing_label    = $args{missing_label} // 'Missing';
    my $csv_file         = $args{csv_file};
    my $csv_missing_label= $args{csv_missing_label} // $missing_label;

    my $rows_ref = print_summary_section(
        title         => $title,
        stats_ref     => $stats_ref,
        missing_label => $missing_label,
    );

    if (defined $csv_file && $csv_file ne '') {
        write_csv_rows($csv_file, $rows_ref, $csv_missing_label, 1); # summary CSV includes TOTAL
    }

    return $rows_ref;
}

sub is_allele_level {
    my ($value) = @_;

    return 0 unless defined $value;

    $value =~ s/^\s+|\s+$//g;
    return 0 if $value eq '';

    # 2-4 colon-separated numeric fields, optional trailing letter at the end
    return $value =~ /^\d+(?::\d+){1,3}[A-Z]?$/;
}

sub debug_serology_donor {
    my ($donor_id, $age, $fields_ref, $col_idx_ref) = @_;
    our (@locus_list, %loci_ser, %loci_dna);

    my @fields  = @$fields_ref;
    my %col_idx = %{$col_idx_ref};

    print "  $donor_id (age $age): serology-only loci\n";

    foreach my $locus (@locus_list) {
        my ($ser_field_1, $ser_field_2) = @{ $loci_ser{$locus} };
        my ($dna_field_1, $dna_field_2) = @{ $loci_dna{$locus} };

        my $ser_val_1 = &normalize_field($fields[ $col_idx{$ser_field_1} ]);
        my $ser_val_2 = &normalize_field($fields[ $col_idx{$ser_field_2} ]);
        my $dna_val_1 = &normalize_field($fields[ $col_idx{$dna_field_1} ]);
        my $dna_val_2 = &normalize_field($fields[ $col_idx{$dna_field_2} ]);

        my $ser_display =
            ($ser_val_1 ne '' || $ser_val_2 ne '')
                ? (($ser_val_1 ne '' ? $ser_val_1 : '-') . ',' . ($ser_val_2 ne '' ? $ser_val_2 : '-'))
                : '-,-';

        my $dna_display =
            ($dna_val_1 ne '' || $dna_val_2 ne '')
                ? (($dna_val_1 ne '' ? $dna_val_1 : '-') . ',' . ($dna_val_2 ne '' ? $dna_val_2 : '-'))
                : '-,-';

        my $has_ser = ($ser_val_1 ne '' || $ser_val_2 ne '');
        my $has_dna = ($dna_val_1 ne '' || $dna_val_2 ne '');

        my $marker = ($has_ser && !$has_dna) ? ' [SERO-ONLY]' : '';
        printf "    %-3s Ser=%-15s DNA=%s%s\n", $locus . ':', $ser_display, $dna_display, $marker;
    }
}

sub debug_allele_donor {
    my ($donor_id, $age, $fields_ref, $col_idx_ref) = @_;
    our (@locus_list, %loci_dna);

    my @fields  = @$fields_ref;
    my %col_idx = %{$col_idx_ref};

    print "  $donor_id (age $age): allele-level resolution lacking\n";

    foreach my $locus (@locus_list) {
        my ($dna_field_1, $dna_field_2) = @{ $loci_dna{$locus} };

        my $dna_val_1 = &normalize_field($fields[$col_idx{$dna_field_1}]);
        my $dna_val_2 = &normalize_field($fields[$col_idx{$dna_field_2}]);

        my $val_1_ok = is_allele_level($dna_val_1);
        my $val_2_ok = ($dna_val_2 eq '') || is_allele_level($dna_val_2);

        my $display_1 = ($dna_val_1 ne '') ? $dna_val_1 : '-';
        my $display_2 = ($dna_val_2 ne '') ? $dna_val_2 : '-';

        my $reason = '';
        if (!$val_1_ok) {
            $reason = " [1st invalid: '$display_1']";
        } elsif (!$val_2_ok) {
            $reason = " [2nd invalid: '$display_2']";
        }
        printf "    %-5s %s, %s%s\n", $locus . ':', $display_1, $display_2, $reason;
    }
}
