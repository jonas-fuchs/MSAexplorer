// This contains some helper functions (written with ChatGPT because JS is a pain)
// to ensure custom behaviour for MSAexplorer
// If anyone can do this better, please consider doing a PR


// get the correct window size for pdf plotting
$(document).on('shiny:connected', function(e) {
    var dimensions = {
        width: $(window).width(),
        height: $(window).height()
    };
    Shiny.setInputValue("window_dimensions_plot", dimensions);
});

// double click for fullscreen plot
document.addEventListener('DOMContentLoaded', function () {
            const plot = document.getElementById('msa_plot');

            // Toggle fullscreen mode on double-click
            plot.ondblclick = function () {
                if (!document.fullscreenElement && !document.webkitFullscreenElement && !document.msFullscreenElement) {
                    // Enter fullscreen
                    if (plot.requestFullscreen) {
                        plot.requestFullscreen();
                    } else if (plot.webkitRequestFullscreen) {
                        plot.webkitRequestFullscreen();
                    } else if (plot.msRequestFullscreen) {
                        plot.msRequestFullscreen();
                    }
                } else {
                    // Exit fullscreen
                    if (document.exitFullscreen) {
                        document.exitFullscreen();
                    } else if (document.webkitExitFullscreen) {
                        document.webkitExitFullscreen();
                    } else if (document.msExitFullscreen) {
                        document.msExitFullscreen();
                    }
                }
            };
        });

// get the correct window size for layout
let resizeTimeout;

$(window).on('resize', function() {
    clearTimeout(resizeTimeout);

    resizeTimeout = setTimeout(function() {
        let dimensions = {
            width: $(window).width(),
            height: $(window).height()
        };
        Shiny.setInputValue("window_dimensions", dimensions);
    }, 500);
});


// Listen for the custom message from the server to update the plot container height
Shiny.addCustomMessageHandler("update-plot-container-height", function(message) {
    // Find the plot container by its ID
    var plotContainer = document.getElementById("msa_plot");
    if (plotContainer) {
        // Update the CSS height property for the plot container
        plotContainer.style.setProperty("height", message.height);
    }
});

// Navbar hides after a while
document.addEventListener("DOMContentLoaded", function () {
    let lastScrollY = window.scrollY;
    const navbar = document.querySelector(".navbar");
    let hideTimeout;
    let isNavbarHidden = false;

    window.addEventListener("scroll", function () {
        clearTimeout(hideTimeout);

        if (window.scrollY > lastScrollY && window.scrollY > 30) {
            if (!isNavbarHidden) {
                hideTimeout = setTimeout(() => {
                    navbar.style.transform = "translateY(-100%)";
                    isNavbarHidden = true;
                }, 500);
            }
        } else if (window.scrollY < lastScrollY && window.scrollY < 100) {
            navbar.style.transform = "translateY(0)";
            isNavbarHidden = false;
        }

        lastScrollY = window.scrollY;
    });
});

// custom sidebar behaviour
function toggleSidebar() {
    var sidebar = document.getElementById("overlay-sidebar");
    var bg = document.getElementById("overlay-bg");
    if (sidebar.classList.contains("show")) {
        sidebar.classList.remove("show");
        bg.style.display = "none";
    } else {
        sidebar.classList.add("show");
        bg.style.display = "block";
    }
}

// collapse navbar on mobile devices
document.addEventListener("click", function(event) {
    var navbar = document.querySelector(".navbar-collapse");
    var toggle = document.querySelector(".navbar-toggler");

    if (!navbar.contains(event.target) && !toggle.contains(event.target)) {
        navbar.classList.remove("show");  // Bootstrap collapse class
    }
});

// Function to update Shiny input when a color changes
function updateColor(inputId, color) {
    if (window.Shiny) {
        Shiny.setInputValue(inputId, color);
    }
}

// Send initial colors to Shiny when the document is ready
document.addEventListener("DOMContentLoaded", function () {
    // Find all color pickers on the page
    const colorPickers = document.querySelectorAll('input[type="color"]');

    // Iterate over each color picker
    colorPickers.forEach(function (picker) {
        const inputId = picker.id; // Use the element's ID as the Shiny input ID
        const initialColor = picker.value;

        // Update Shiny with the initial value
        const checkShinyReady = setInterval(function () {
            if (window.Shiny && typeof Shiny.setInputValue === 'function') {
                clearInterval(checkShinyReady);
                updateColor(inputId, initialColor);
            }
        }, 100);

        // Add event listener for changes
        picker.addEventListener('change', function () {
            updateColor(inputId, picker.value);
        });
    });
});