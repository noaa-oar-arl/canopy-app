// NOAA Canopy-App Documentation Enhancements

document.addEventListener('DOMContentLoaded', function() {
    // Add NOAA branding enhancements
    addNOAABranding();

    // Add scroll animations
    addScrollAnimations();

    // Add enhanced search functionality
    enhanceSearch();

    // Add copy code button enhancements
    enhanceCodeBlocks();
});

function addNOAABranding() {
    // Add NOAA badge to footer
    const footer = document.querySelector('.md-footer-meta');
    if (footer) {
        const noaaBadge = document.createElement('div');
        noaaBadge.innerHTML = `
            <div style="text-align: center; margin-top: 1rem; padding: 1rem; border-top: 1px solid rgba(255,255,255,0.1);">
                <p style="margin: 0; opacity: 0.8; font-size: 0.9rem;">
                    🌊 Developed by <strong>NOAA Air Resources Laboratory</strong> 🌍
                </p>
                <p style="margin: 0.5rem 0 0 0; opacity: 0.6; font-size: 0.8rem;">
                    Supporting atmospheric research and environmental protection
                </p>
            </div>
        `;
        footer.appendChild(noaaBadge);
    }

    // Add atmospheric icons to navigation if on homepage
    if (window.location.pathname.endsWith('/') || window.location.pathname.includes('index')) {
        addAtmosphericIcons();
    }
}

function addAtmosphericIcons() {
    // Add floating atmospheric elements for visual appeal
    const content = document.querySelector('.md-content');
    if (content) {
        const atmosphericBg = document.createElement('div');
        atmosphericBg.className = 'atmospheric-background';
        atmosphericBg.innerHTML = `
            <div class="cloud cloud1">☁️</div>
            <div class="cloud cloud2">🌤️</div>
            <div class="cloud cloud3">☁️</div>
            <div class="leaf leaf1">🍃</div>
            <div class="leaf leaf2">🌿</div>
        `;
        content.appendChild(atmosphericBg);

        // Add CSS for floating elements
        const style = document.createElement('style');
        style.textContent = `
            .atmospheric-background {
                position: fixed;
                top: 0;
                left: 0;
                width: 100%;
                height: 100%;
                pointer-events: none;
                z-index: 0;
                opacity: 0.1;
            }

            .cloud, .leaf {
                position: absolute;
                font-size: 2rem;
                animation: float 20s infinite linear;
            }

            .cloud1 { top: 10%; left: 10%; animation-delay: 0s; }
            .cloud2 { top: 20%; right: 15%; animation-delay: -5s; }
            .cloud3 { top: 30%; left: 70%; animation-delay: -10s; }
            .leaf1 { bottom: 20%; left: 20%; animation-delay: -8s; animation-duration: 15s; }
            .leaf2 { bottom: 30%; right: 25%; animation-delay: -12s; animation-duration: 18s; }

            @keyframes float {
                0% { transform: translateY(0px) rotate(0deg); }
                50% { transform: translateY(-20px) rotate(180deg); }
                100% { transform: translateY(0px) rotate(360deg); }
            }

            @media (max-width: 768px) {
                .atmospheric-background { display: none; }
            }
        `;
        document.head.appendChild(style);
    }
}

function addScrollAnimations() {
    // Add scroll-triggered animations for better user experience
    const observerOptions = {
        threshold: 0.1,
        rootMargin: '0px 0px -50px 0px'
    };

    const observer = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                entry.target.style.opacity = '1';
                entry.target.style.transform = 'translateY(0)';
            }
        });
    }, observerOptions);

    // Animate elements as they come into view
    document.querySelectorAll('.md-typeset h2, .md-typeset h3, .admonition, table').forEach(el => {
        el.style.opacity = '0';
        el.style.transform = 'translateY(20px)';
        el.style.transition = 'opacity 0.6s ease, transform 0.6s ease';
        observer.observe(el);
    });
}

function enhanceSearch() {
    // Add weather-related search suggestions
    const searchInput = document.querySelector('.md-search__input');
    if (searchInput) {
        const suggestions = [
            'wind profile', 'biogenic emissions', 'dry deposition',
            'photolysis', 'canopy height', 'LAI', 'meteorology',
            'MEGAN', 'namelist', 'configuration', 'installation'
        ];

        searchInput.addEventListener('focus', () => {
            if (!searchInput.value) {
                const randomSuggestion = suggestions[Math.floor(Math.random() * suggestions.length)];
                searchInput.placeholder = `Try searching for "${randomSuggestion}"...`;
            }
        });

        searchInput.addEventListener('blur', () => {
            searchInput.placeholder = 'Search documentation...';
        });
    }
}

function enhanceCodeBlocks() {
    // Add enhanced copy functionality for code blocks
    document.querySelectorAll('.highlight').forEach(block => {
        const copyButton = block.querySelector('.md-clipboard');
        if (copyButton) {
            copyButton.addEventListener('click', () => {
                // Add visual feedback
                const originalText = copyButton.textContent;
                copyButton.textContent = '✅ Copied!';
                copyButton.style.color = '#4CAF50';

                setTimeout(() => {
                    copyButton.textContent = originalText;
                    copyButton.style.color = '';
                }, 2000);
            });
        }
    });
}

// Add keyboard shortcuts for better navigation
document.addEventListener('keydown', (e) => {
    // Alt + H for homepage
    if (e.altKey && e.key === 'h') {
        e.preventDefault();
        window.location.href = '/';
    }

    // Alt + S for search
    if (e.altKey && e.key === 's') {
        e.preventDefault();
        const searchInput = document.querySelector('.md-search__input');
        if (searchInput) {
            searchInput.focus();
        }
    }
});

// Add theme transition effects
const themeToggle = document.querySelector('[data-md-component="palette"]');
if (themeToggle) {
    themeToggle.addEventListener('change', () => {
        document.body.style.transition = 'all 0.3s ease';
        setTimeout(() => {
            document.body.style.transition = '';
        }, 300);
    });
}
