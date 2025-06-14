#!/bin/bash
# Comprehensive documentation management script for Canopy-App
# Handles setup, building, serving, and deployment
# Usage: ./scripts/docs.sh [command] [options]

set -e

# Configuration
DOCS_DIR="docs"
SITE_DIR="site"
REQUIREMENTS_FILE="requirements-docs.txt"
DEFAULT_VERSION="latest"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Functions
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
    echo "║              Canopy-App Documentation Manager               ║"
    echo "║         Setup • Build • Serve • Deploy • Manage             ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

show_help() {
    echo "Canopy-App Documentation Manager"
    echo ""
    echo "USAGE:"
    echo "  $0 <command> [options]"
    echo ""
    echo "COMMANDS:"
    echo "  setup              - First-time setup and validation"
    echo "  build [version]    - Build documentation"
    echo "  serve [version]    - Serve documentation locally"
    echo "  deploy [version]   - Deploy to GitHub Pages"
    echo "  list               - List available versions"
    echo "  clean              - Clean build artifacts"
    echo "  help               - Show this help message"
    echo ""
    echo "EXAMPLES:"
    echo "  $0 setup                    # First-time setup"
    echo "  $0 serve                    # Serve locally with hot reload"
    echo "  $0 build                    # Build documentation"
    echo "  $0 deploy                   # Deploy latest version"
    echo "  $0 deploy v1.0.0            # Deploy specific version"
    echo ""
}

# Check if we're in the right directory
check_directory() {
    if [ ! -f "mkdocs.yml" ]; then
        log_error "mkdocs.yml not found. Please run this script from the project root."
        exit 1
    fi
}

# Check requirements and dependencies
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

# Install documentation dependencies
install_dependencies() {
    log_info "Installing documentation dependencies..."
    
    if [ -f "$REQUIREMENTS_FILE" ]; then
        pip install -r "$REQUIREMENTS_FILE"
        log_success "Documentation dependencies installed"
    else
        log_error "requirements-docs.txt not found"
        exit 1
    fi
}

# Check GitHub Pages configuration
check_github_pages() {
    log_info "Checking GitHub Pages configuration..."
    
    # Check if GitHub CLI is available
    if command -v gh &> /dev/null; then
        log_info "GitHub CLI detected, checking repository settings..."
        
        # Try to get repository info
        if gh repo view noaa-oar-arl/canopy-app &> /dev/null; then
            log_success "Repository accessible via GitHub CLI"
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

# Build documentation
build_docs() {
    log_info "Building documentation..."
    
    # Clean previous build
    if [ -d "$SITE_DIR" ]; then
        rm -rf "$SITE_DIR"
        log_info "Cleaned previous build"
    fi
    
    # Build with MkDocs
    mkdocs build --verbose --clean
    log_success "Documentation built successfully"
    
    # Show build stats
    if [ -d "$SITE_DIR" ]; then
        SITE_SIZE=$(du -sh "$SITE_DIR" | cut -f1)
        PAGE_COUNT=$(find "$SITE_DIR" -name "*.html" | wc -l)
        log_info "Generated site size: $SITE_SIZE"
        log_info "Generated pages: $PAGE_COUNT"
    fi
}

# Serve documentation locally
serve_docs() {
    log_info "Starting local documentation server..."
    log_info "Documentation will be available at: http://127.0.0.1:8000"
    log_info "Press Ctrl+C to stop the server"
    
    mkdocs serve
}

# Deploy to GitHub Pages
deploy_docs() {
    local VERSION=${1:-$DEFAULT_VERSION}
    log_info "Deploying documentation to GitHub Pages (version: $VERSION)..."
    
    # Check if git is configured
    if ! git config user.name &> /dev/null || ! git config user.email &> /dev/null; then
        log_warning "Git user not configured. Setting default values..."
        git config --local user.email "docs@noaa-oar-arl.github.io"
        git config --local user.name "Documentation Bot"
    fi
    
    # Deploy with mike for versioning
    if command -v mike &> /dev/null; then
        log_info "Deploying version '$VERSION' with mike..."
        mike deploy --push --update-aliases "$VERSION"
        
        if [ "$VERSION" == "latest" ] || [ "$VERSION" == "main" ]; then
            mike set-default --push "$VERSION"
            log_info "Set '$VERSION' as default version"
        fi
        
        log_success "Documentation deployed successfully with versioning"
        log_info "Available at: https://noaa-oar-arl.github.io/canopy-app/"
    else
        log_warning "mike not found. Deploying without versioning..."
        mkdocs gh-deploy --force
        log_success "Documentation deployed successfully"
        log_info "Available at: https://noaa-oar-arl.github.io/canopy-app/"
    fi
}

# List available documentation versions
list_versions() {
    log_info "Available documentation versions:"
    
    if command -v mike &> /dev/null; then
        mike list
    else
        log_warning "mike not installed. Cannot list versions."
        log_info "Install mike for version management: pip install mike"
    fi
}

# Clean build artifacts
clean_docs() {
    log_info "Cleaning documentation build artifacts..."
    
    if [ -d "$SITE_DIR" ]; then
        rm -rf "$SITE_DIR"
        log_success "Removed $SITE_DIR directory"
    fi
    
    # Clean any Python cache
    find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
    find . -name "*.pyc" -delete 2>/dev/null || true
    
    log_success "Cleanup completed"
}

# Test documentation build
test_docs() {
    log_info "Testing documentation build..."
    
    # Test MkDocs build
    if mkdocs build --clean --quiet; then
        log_success "Documentation builds successfully"
        
        # Check if site directory was created
        if [ -d "$SITE_DIR" ]; then
            SITE_SIZE=$(du -sh "$SITE_DIR" | cut -f1)
            log_info "Generated site size: $SITE_SIZE"
            
            # Count pages
            PAGE_COUNT=$(find "$SITE_DIR" -name "*.html" | wc -l)
            log_info "Generated pages: $PAGE_COUNT"
        fi
    else
        log_error "Documentation build failed"
        log_error "Check mkdocs.yml configuration and run 'mkdocs build' for details"
        exit 1
    fi
}

# Setup command - first-time setup and validation
setup_command() {
    print_banner
    check_requirements
    install_dependencies
    check_github_pages
    test_docs
    
    echo
    log_success "Documentation setup completed successfully!"
    echo
    echo -e "${BLUE}Next Steps:${NC}"
    echo "  1. 📝 Edit documentation files in the 'docs/' directory"
    echo "  2. 🔧 Customize 'mkdocs.yml' configuration if needed"
    echo "  3. 🧪 Test locally: $0 serve"
    echo "  4. 📤 Deploy: $0 deploy"
    echo "  5. 🌐 Visit: https://noaa-oar-arl.github.io/canopy-app/"
    echo
    echo -e "${YELLOW}Useful Commands:${NC}"
    echo "  • Local development: $0 serve"
    echo "  • Build only: $0 build"
    echo "  • Deploy: $0 deploy"
    echo "  • Clean: $0 clean"
    echo
}

# Main execution
main() {
    check_directory
    
    local COMMAND=${1:-"help"}
    local VERSION=${2:-$DEFAULT_VERSION}
    
    case $COMMAND in
        "setup")
            setup_command
            ;;
        "build")
            check_requirements
            install_dependencies
            build_docs
            ;;
        "serve")
            check_requirements
            install_dependencies
            serve_docs
            ;;
        "deploy")
            check_requirements
            install_dependencies
            build_docs
            deploy_docs "$VERSION"
            ;;
        "list")
            list_versions
            ;;
        "clean")
            clean_docs
            ;;
        "test")
            check_requirements
            install_dependencies
            test_docs
            ;;
        "help"|"-h"|"--help")
            show_help
            ;;
        *)
            log_error "Unknown command: $COMMAND"
            echo
            show_help
            exit 1
            ;;
    esac
}

main "$@"
