# Setup Instructions for GitHub

This document explains how to set up this repository on GitHub.

## Initial Setup

### 1. Create a new repository on GitHub

1. Go to [GitHub](https://github.com) and log in
2. Click the "+" icon in the top right and select "New repository"
3. Choose a repository name (e.g., `gcd-analyzer` or `hla-typing-analyzer`)
4. Add a description: "HLA typing resolution quality analysis tool for ION CSV data"
5. Choose "Public" visibility
6. **Do NOT** initialize with README, .gitignore, or LICENSE (we already have these)
7. Click "Create repository"

### 2. Initialize and push from command line

```bash
# Navigate to the repository directory
cd /path/to/gcd-analyzer

# Initialize git repository
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: GCD Analyzer v1.0"

# Add the remote repository (replace <username> with your GitHub username)
git remote add origin https://github.com/<username>/gcd-analyzer.git

# Push to GitHub
git branch -M main
git push -u origin main
```

### 3. Verify on GitHub

After pushing, verify that:
- All files are visible on GitHub
- README.md is displayed on the repository home page
- Examples are properly organized

## Repository Structure

The repository should contain:

```
gcd-analyzer/
├── analyze_gcd.pl           # Main analysis script
├── README.md                # Documentation
├── LICENSE                  # MIT License
├── .gitignore              # Git ignore rules
├── SETUP.md                # This file
└── examples/
    ├── input/              # Example input files
    │   ├── ION-8888-D.csv
    │   └── ION-9999-D.csv
    └── output/             # Example output files
        ├── ION-8888-D.csv
        ├── ION-9999-D.csv
        ├── summary_C-Donors.csv
        ├── summary_D-Donors.csv
        └── summary_Overall.csv
```

## Adding Topics and Tags

To make your repository more discoverable on GitHub:

1. Go to your repository on GitHub
2. Click the gear icon (⚙️) next to "About"
3. Add topics such as:
   - `hla-typing`
   - `bioinformatics`
   - `perl`
   - `data-analysis`
   - `medical-research`
   - `transplantation`

## Optional: GitHub Actions

You could add a simple GitHub Action to test the script. Create `.github/workflows/test.yml`:

```yaml
name: Test

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Install Perl
        run: sudo apt-get install -y perl
      - name: Test script
        run: |
          chmod +x analyze_gcd.pl
          ./analyze_gcd.pl --dir examples/input
          ./analyze_gcd.pl --help
```

## Sharing with Colleagues

Once the repository is set up:

1. Share the repository URL: `https://github.com/<username>/gcd-analyzer`
2. Colleagues can clone it with: `git clone https://github.com/<username>/gcd-analyzer.git`
3. They can report issues via the GitHub Issues tab
4. They can contribute via Pull Requests

## Updating the Repository

To update after making changes:

```bash
# Check status
git status

# Add changes
git add .

# Commit changes
git commit -m "Description of changes"

# Push to GitHub
git push
```

## Additional Resources

- [GitHub Documentation](https://docs.github.com/)
- [Git Tutorial](https://git-scm.com/docs/gittutorial)
- [Markdown Guide](https://www.markdownguide.org/)
