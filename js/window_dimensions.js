// includes short js script to get the correct window size for pdf plotting

$(document).on('shiny:connected', function(e) {
    var dimensions = {
        width: $(window).width(),
        height: $(window).height()
    };
    Shiny.setInputValue("window_dimensions", dimensions);
});

$(window).on('resize', function(e) {
    var dimensions = {
        width: $(window).width(),
        height: $(window).height()
    };
    Shiny.setInputValue("window_dimensions", dimensions);
});