#!/bin/bash
# Setup script for Canopy-App GitHub Pages documentation
# This script helps configure the repository for automatic documentation deployment

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_banner() {
    echo -e "${BLUE}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║              Canopy-App Documentation Setup                 ║"
    echo "║         GitHub Pages + MkDocs Configuration Tool            ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

check_requirements() {
    log_info "Checking requirements..."

    # Check if we're in a git repository
    if ! git rev-parse --git-dir > /dev/null 2>&1; then
        log_error "This script must be run from within a git repository"
        exit 1
    fi

    # Check if we're in the right repository
    REPO_URL=$(git config --get remote.origin.url 2>/dev/null || echo "")
    if [[ ! "$REPO_URL" =~ "noaa-oar-arl/canopy-app" ]]; then
        log_warning "This script is designed for the noaa-oar-arl/canopy-app repository"
        log_warning "Current repository: $REPO_URL"
        read -p "Continue anyway? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            exit 1
        fi
    fi

    # Check if Python is available
    if ! command -v python3 &> /dev/null; then
        log_error "Python 3 is required but not installed"
        exit 1
    fi

    # Check if pip is available
    if ! command -v pip &> /dev/null; then
        log_error "pip is required but not installed"
        exit 1
    fi

    log_success "Requirements check passed"
}

install_dependencies() {
    log_info "Installing documentation dependencies..."

    if [ -f "requirements-docs.txt" ]; then
        pip install -r requirements-docs.txt
        log_success "Documentation dependencies installed"
    else
        log_error "requirements-docs.txt not found"
        exit 1
    fi
}

check_github_pages() {
    log_info "Checking GitHub Pages configuration..."

    # Check if GitHub CLI is available
    if command -v gh &> /dev/null; then
        log_info "GitHub CLI detected, checking repository settings..."

        # Try to get repository info
        if gh repo view noaa-oar-arl/canopy-app &> /dev/null; then
            log_success "Repository accessible via GitHub CLI"

            # Check if Pages is enabled (this might require admin access)
            log_info "Note: GitHub Pages should be configured to deploy from GitHub Actions"
            log_info "Visit: https://github.com/noaa-oar-arl/canopy-app/settings/pages"
        else
            log_warning "Cannot access repository via GitHub CLI"
            log_warning "You may need to authenticate: gh auth login"
        fi
    else
        log_warning "GitHub CLI not installed"
        log_info "Consider installing it for easier repository management: https://cli.github.com/"
    fi

    log_info "Manual GitHub Pages setup:"
    echo "  1. Go to: https://github.com/noaa-oar-arl/canopy-app/settings/pages"
    echo "  2. Set Source to 'GitHub Actions'"
    echo "  3. Save the configuration"
}

test_documentation() {
    log_info "Testing documentation build..."

    # Test MkDocs build
    if mkdocs build --clean --quiet; then
        log_success "Documentation builds successfully"

        # Check if site directory was created
        if [ -d "site" ]; then
            SITE_SIZE=$(du -sh site | cut -f1)
            log_info "Generated site size: $SITE_SIZE"

            # Count pages
            PAGE_COUNT=$(find site -name "*.html" | wc -l)
            log_info "Generated pages: $PAGE_COUNT"
        fi
    else
        log_error "Documentation build failed"
        log_error "Check mkdocs.yml configuration and run 'mkdocs build' for details"
        exit 1
    fi
}

show_next_steps() {
    echo
    log_success "Setup completed successfully!"
    echo
    echo -e "${BLUE}Next Steps:${NC}"
    echo "  1. 📝 Edit documentation files in the 'docs/' directory"
    echo "  2. 🔧 Customize 'mkdocs.yml' configuration if needed"
    echo "  3. 🧪 Test locally with: ./scripts/deploy-docs.sh latest serve"
    echo "  4. 📤 Push changes to trigger automatic deployment"
    echo "  5. 🌐 Visit your documentation at: https://noaa-oar-arl.github.io/canopy-app/"
    echo
    echo -e "${YELLOW}Useful Commands:${NC}"
    echo "  • Local development: mkdocs serve"
    echo "  • Build documentation: mkdocs build"
    echo "  • Deploy manually: ./scripts/deploy-docs.sh"
    echo "  • Check workflows: https://github.com/noaa-oar-arl/canopy-app/actions"
    echo
}

main() {
    print_banner

    check_requirements
    install_dependencies
    check_github_pages
    test_documentation
    show_next_steps
}

# Run main function
main "$@"
