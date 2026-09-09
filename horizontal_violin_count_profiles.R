# horizontal count-scaled profiles
# Compatible with R 4.1.1; no additional packages are required.
#
# IMPORTANT ENCODING
# At every observed repeat, the FULL vertical thickness of the profile equals
# the neuronal count in SVG pixels. Example: 40 neurons = 40 SVG pixels.
# This is therefore a count-scaled violin-shaped profile, not a KDE violin.


# =============================================================================
# 1. SETTINGS
# =============================================================================

GROUP_ORDER <- c("Wt", "Homo", "Het")

COLORS <- c(
 Wt = "#7A7A7A",
 Homo = "#E67E22",
 Het = "#55ABBF"
)

TITLE <- "TH in SNpc"
SCRIPT_VERSION <- "compact_v3"

# Output folder. Defaults to a local "output" directory (created if missing).
# Override at runtime, e.g.:
#   Rscript horizontal_violin_count_profiles.R --output /path/to/figures
OUTPUT_DIR <- file.path(getwd(), "output")
args <- commandArgs(trailingOnly = TRUE)
output_flag <- match("--output", args)
if (!is.na(output_flag) && output_flag < length(args)) {
  OUTPUT_DIR <- args[output_flag + 1]
}

OUTPUT_SVG <- file.path(
 OUTPUT_DIR,
 "TH_SNc.svg"
)

# Exact scale. Do not change this if 1 neuron must equal 1 SVG pixel.
PIXELS_PER_NEURON <- 1

# Horizontal geometry.
MAX_REPEAT <- 11
REPEAT_GAP_PX <- 40
TAIL_LENGTH_PX <- 18
END_THICKNESS_PX <- 4
CAP_LENGTH_PX <- 5

# Center-to-center vertical distance between adjacent animals in a group.
# Smaller values put the violins closer. Try 60-90 px for this dataset.
# The same value is used for Wt, Homo and Het.
ANIMAL_CENTER_GAP_PX <- 150

# Appearance.
VIOLIN_FILL_OPACITY <- 0.36
VIOLIN_STROKE_OPACITY <- 0.92
VIOLIN_STROKE_WIDTH_PX <- 1.0
POINT_RADIUS_PX <- 2.4
POINT_OPACITY <- 0.78

# Set FALSE if you want only points and profiles, without printed counts.
SHOW_COUNT_LABELS <- TRUE

# Opens the finished SVG in your default browser after it is saved.
OPEN_AFTER_SAVE <- TRUE


# =============================================================================
# 2. DATA
# =============================================================================
# Each inner vector holds the TH+ neuron counts across the 11 coronal
# sections (repeat positions) for one animal.

DATA <- list(
 Wt = list(
 c(63, 167, 168, 125, 96, 106, 121, 126, 88, 88, 103),
 c(78, 157, 146, 159, 184, 166, 145, 110, 103, 138, 136),
 c(67, 119, 142, 182, 170, 177, 120, 100, 120, 70, 84),
 c(27, 88, 196, 169, 156, 136, 96, 123, 113, 104, 46),
 c(75, 118, 184, 197, 153, 184, 182, 177, 131, 120, 74),
 c(136, 125, 146, 185, 165, 135, 125, 120, 135, 126, 99)
 ),

 Homo = list(
 c(23, 68, 84, 88, 100, 95, 77, 77, 90, 100, 57),
 c(60, 83, 115, 118, 84, 118, 75, 89, 108, 98, 58),
 c(60, 80, 106, 135, 121, 65, 123, 54, 82, 84, 119),
 c(50, 58, 93, 110, 91, 90, 108, 91, 78, 57, 91),
 c(63, 79, 111, 140, 127, 129, 65, 89, 83, 75, 115)
 ),

 Het = list(
 c(81, 113, 117, 111, 103, 101, 94, 113, 61, 67, 60),
 c(54, 71, 96, 136, 150, 123, 116, 87, 39, 60, 99),
 c(54, 65, 115, 147, 139, 105, 115, 96, 125, 111, 140),
 c(65, 107, 134, 161, 175, 119, 117, 119, 152, 121, 128),
 c(35, 75, 97, 149, 132, 151, 117, 108, 76, 104, 116),
 c(100, 127, 143, 142, 132, 136, 91, 78, 84, 136, 109),
 c(101, 104, 149, 156, 125, 139, 115, 107, 85, 107, 124),
 c(70, 108, 100, 104, 75, 82, 67, 89, 73, 81, 97),
 c(22, 24, 32, 45, 35, 21, 24, 29, 16, 57, 23),
 c(70, 107, 153, 153, 160, 127, 128, 121, 92, 147, 150)
 )
)

# Fail immediately if animal identifiers are accidentally reintroduced.
stopifnot(all(vapply(DATA, function(group) is.null(names(group)), logical(1))))


# =============================================================================
# 3. SVG AND SMOOTH-CURVE HELPERS
# =============================================================================

xml_escape <- function(text) {
 text <- gsub("&", "&amp;", text, fixed = TRUE)
 text <- gsub("<", "&lt;", text, fixed = TRUE)
 text <- gsub(">", "&gt;", text, fixed = TRUE)
 text <- gsub('"', "&quot;", text, fixed = TRUE)
 text
}


number <- function(x) {
 formatC(x, digits = 3, format = "f")
}


# Shape-preserving PCHIP derivatives. These produce smooth cubic Bezier curves
# that pass through every measured thickness without sharp polygonal corners.
pchip_slopes <- function(x, y) {
 n <- length(x)

 if (n == 2) {
 slope <- (y[2] - y[1]) / (x[2] - x[1])
 return(c(slope, slope))
 }

 h <- diff(x)
 delta <- diff(y) / h
 derivative <- numeric(n)

 # Shape-preserving endpoint derivatives.
 derivative[1] <-
 ((2 * h[1] + h[2]) * delta[1] - h[1] * delta[2]) /
 (h[1] + h[2])

 if (sign(derivative[1]) != sign(delta[1])) {
 derivative[1] <- 0
 } else if (
 sign(delta[1]) != sign(delta[2]) &&
 abs(derivative[1]) > abs(3 * delta[1])
 ) {
 derivative[1] <- 3 * delta[1]
 }

 derivative[n] <-
 ((2 * h[n - 1] + h[n - 2]) * delta[n - 1] -
 h[n - 1] * delta[n - 2]) /
 (h[n - 1] + h[n - 2])

 if (sign(derivative[n]) != sign(delta[n - 1])) {
 derivative[n] <- 0
 } else if (
 sign(delta[n - 1]) != sign(delta[n - 2]) &&
 abs(derivative[n]) > abs(3 * delta[n - 1])
 ) {
 derivative[n] <- 3 * delta[n - 1]
 }

 # Weighted harmonic means prevent overshoot at internal anchors.
 for (index in 2:(n - 1)) {
 previous_delta <- delta[index - 1]
 next_delta <- delta[index]

 if (
 previous_delta == 0 ||
 next_delta == 0 ||
 sign(previous_delta) != sign(next_delta)
 ) {
 derivative[index] <- 0
 } else {
 weight_1 <- 2 * h[index] + h[index - 1]
 weight_2 <- h[index] + 2 * h[index - 1]

 derivative[index] <-
 (weight_1 + weight_2) /
 (
 weight_1 / previous_delta +
 weight_2 / next_delta
 )
 }
 }

 derivative
}


bezier_segments <- function(x, y, derivative) {
 segments <- character(length(x) - 1)

 for (index in seq_len(length(x) - 1)) {
 interval <- x[index + 1] - x[index]

 control_1_x <- x[index] + interval / 3
 control_1_y <- y[index] + derivative[index] * interval / 3

 control_2_x <- x[index + 1] - interval / 3
 control_2_y <-
 y[index + 1] - derivative[index + 1] * interval / 3

 segments[index] <- paste(
 "C",
 number(control_1_x), number(control_1_y),
 number(control_2_x), number(control_2_y),
 number(x[index + 1]), number(y[index + 1])
 )
 }

 segments
}


# Produces one true cubic-Bezier SVG path. At an observed repeat position,
# bottom_y - top_y equals count * PIXELS_PER_NEURON exactly.
violin_path <- function(values, repeat_x, center_y) {
 observed <- which(!is.na(values))

 observed_x <- repeat_x[observed]
 observed_thickness <- values[observed] * PIXELS_PER_NEURON

 x <- c(
 min(observed_x) - TAIL_LENGTH_PX,
 observed_x,
 max(observed_x) + TAIL_LENGTH_PX
 )

 thickness <- c(
 END_THICKNESS_PX,
 observed_thickness,
 END_THICKNESS_PX
 )

 thickness_derivative <- pchip_slopes(x, thickness)

 top_y <- center_y - thickness / 2
 bottom_y <- center_y + thickness / 2

 top_derivative <- -thickness_derivative / 2
 bottom_derivative <- thickness_derivative / 2

 top_segments <- bezier_segments(x, top_y, top_derivative)
 bottom_segments_forward <- bezier_segments(
 x,
 bottom_y,
 bottom_derivative
 )

 # Reverse each lower Bezier segment so the closed outline runs right-to-left.
 bottom_segments_reverse <- character(length(x) - 1)

 for (index in (length(x) - 1):1) {
 interval <- x[index + 1] - x[index]

 forward_control_1_x <- x[index] + interval / 3
 forward_control_1_y <-
 bottom_y[index] + bottom_derivative[index] * interval / 3

 forward_control_2_x <- x[index + 1] - interval / 3
 forward_control_2_y <-
 bottom_y[index + 1] -
 bottom_derivative[index + 1] * interval / 3

 reverse_position <- length(x) - index

 bottom_segments_reverse[reverse_position] <- paste(
 "C",
 number(forward_control_2_x), number(forward_control_2_y),
 number(forward_control_1_x), number(forward_control_1_y),
 number(x[index]), number(bottom_y[index])
 )
 }

 left_x <- x[1]
 right_x <- x[length(x)]
 left_radius <- thickness[1] / 2
 right_radius <- thickness[length(thickness)] / 2
 kappa <- 0.55228475

 # Two cubic curves make each rounded half-ellipse cap.
 right_cap <- c(
 paste(
 "C",
 number(right_x + kappa * CAP_LENGTH_PX),
 number(center_y - right_radius),
 number(right_x + CAP_LENGTH_PX),
 number(center_y - kappa * right_radius),
 number(right_x + CAP_LENGTH_PX),
 number(center_y)
 ),
 paste(
 "C",
 number(right_x + CAP_LENGTH_PX),
 number(center_y + kappa * right_radius),
 number(right_x + kappa * CAP_LENGTH_PX),
 number(center_y + right_radius),
 number(right_x),
 number(center_y + right_radius)
 )
 )

 left_cap <- c(
 paste(
 "C",
 number(left_x - kappa * CAP_LENGTH_PX),
 number(center_y + left_radius),
 number(left_x - CAP_LENGTH_PX),
 number(center_y + kappa * left_radius),
 number(left_x - CAP_LENGTH_PX),
 number(center_y)
 ),
 paste(
 "C",
 number(left_x - CAP_LENGTH_PX),
 number(center_y - kappa * left_radius),
 number(left_x - kappa * CAP_LENGTH_PX),
 number(center_y - left_radius),
 number(left_x),
 number(center_y - left_radius)
 )
 )

 paste(
 "M", number(x[1]), number(top_y[1]),
 paste(top_segments, collapse = " "),
 paste(right_cap, collapse = " "),
 paste(bottom_segments_reverse, collapse = " "),
 paste(left_cap, collapse = " "),
 "Z"
 )
}


# =============================================================================
# 4. CREATE THE EDITABLE SVG
# =============================================================================

create_horizontal_count_plot <- function(output_file = OUTPUT_SVG) {
 if (!dir.exists(dirname(output_file))) {
 dir.create(
 dirname(output_file),
 recursive = TRUE,
 showWarnings = FALSE
 )
 }

 if (!dir.exists(dirname(output_file))) {
 stop("Could not create the folder: ", dirname(output_file))
 }

 all_counts <- unlist(DATA, recursive = TRUE, use.names = FALSE)
 max_count <- max(all_counts, na.rm = TRUE)
 max_animals <- max(vapply(DATA, length, integer(1)))

 animal_center_gap <- ANIMAL_CENTER_GAP_PX

 panel_width <- 500
 panel_gap <- 40
 left_margin <- 35
 right_margin <- 130
 profile_left_in_panel <- 125

 top_margin <- 165
 bottom_margin <- 105

 profile_block_height <-
 (max_animals - 1) * animal_center_gap

 canvas_width <-
 left_margin +
 length(GROUP_ORDER) * panel_width +
 (length(GROUP_ORDER) - 1) * panel_gap +
 right_margin

 canvas_height <-
 top_margin +
 max_count * PIXELS_PER_NEURON +
 profile_block_height +
 bottom_margin

 plot_mid_y <-
 top_margin +
 max_count * PIXELS_PER_NEURON / 2 +
 profile_block_height / 2

 svg <- character()

 add <- function(...) {
 svg <<- c(svg, paste0(...))
 }

 add('<?xml version="1.0" encoding="UTF-8"?>')
 add(
 '<svg xmlns="http://www.w3.org/2000/svg" ',
 'width="', canvas_width, 'px" ',
 'height="', canvas_height, 'px" ',
 'viewBox="0 0 ', canvas_width, ' ', canvas_height, '">'
 )

 add(
 '<!-- ', SCRIPT_VERSION,
 '; animal-center-gap=', ANIMAL_CENTER_GAP_PX,
 'px; sample-labels=disabled -->'
 )

 add(
 '<rect x="0" y="0" width="', canvas_width,
 '" height="', canvas_height, '" fill="#FFFFFF"/>'
 )

 add(
 '<text x="', left_margin,
 '" y="42" font-family="Arial, sans-serif" ',
 'font-size="22" font-weight="700" fill="#171A1F">',
 xml_escape(TITLE),
 '</text>'
 )

 add(
 '<text x="', left_margin,
 '" y="66" font-family="Arial, sans-serif" ',
 'font-size="12" fill="#59616B">',
 'Horizontal count-scaled profiles: full thickness = neuronal count ',
 '(1 neuron = 1 SVG pixel)',
 '</text>'
 )

 # Legend.
 legend_start_x <- canvas_width - 285

 for (group_index in seq_along(GROUP_ORDER)) {
 group <- GROUP_ORDER[group_index]
 legend_x <- legend_start_x + (group_index - 1) * 92

 add(
 '<rect x="', legend_x,
 '" y="30" width="18" height="11" ',
 'fill="', COLORS[[group]], '" fill-opacity="0.45" ',
 'stroke="', COLORS[[group]], '" stroke-width="1"/>'
 )

 add(
 '<text x="', legend_x + 25,
 '" y="40" font-family="Arial, sans-serif" ',
 'font-size="12" fill="#20242A">', group, '</text>'
 )
 }

 for (group_index in seq_along(GROUP_ORDER)) {
 group <- GROUP_ORDER[group_index]
 animals <- DATA[[group]]
 color <- COLORS[[group]]

 panel_left <-
 left_margin + (group_index - 1) * (panel_width + panel_gap)

 repeat_x <-
 panel_left +
 profile_left_in_panel +
 (seq_len(MAX_REPEAT) - 1) * REPEAT_GAP_PX

 group_span <- (length(animals) - 1) * animal_center_gap
 first_animal_y <- plot_mid_y - group_span / 2

 add(
 '<text x="', panel_left + panel_width / 2,
 '" y="105" text-anchor="middle" ',
 'font-family="Arial, sans-serif" font-size="17" ',
 'font-weight="700" fill="', color, '">',
 group,
 '</text>'
 )

 # Repeat grid: the x position is the same for every animal and group.
 plot_top <- top_margin - 20
 plot_bottom <- canvas_height - bottom_margin + 15

 for (repeat_index in seq_len(MAX_REPEAT)) {
 add(
 '<line x1="', repeat_x[repeat_index],
 '" y1="', plot_top,
 '" x2="', repeat_x[repeat_index],
 '" y2="', plot_bottom,
 '" stroke="#E2E6EB" stroke-width="0.75"/>'
 )

 add(
 '<text x="', repeat_x[repeat_index],
 '" y="', canvas_height - 66,
 '" text-anchor="middle" font-family="Arial, sans-serif" ',
 'font-size="11" fill="#59616B">',
 repeat_index,
 '</text>'
 )
 }

 add(
 '<text x="', mean(range(repeat_x)),
 '" y="', canvas_height - 37,
 '" text-anchor="middle" font-family="Arial, sans-serif" ',
 'font-size="12" fill="#20242A">Coronal section</text>'
 )

 for (animal_index in seq_along(animals)) {
 values <- animals[[animal_index]]
 center_y <- first_animal_y + (animal_index - 1) * animal_center_gap

 path_data <- violin_path(values, repeat_x, center_y)

 add(
 '<path d="', path_data,
 '" fill="', color,
 '" fill-opacity="', VIOLIN_FILL_OPACITY,
 '" stroke="', color,
 '" stroke-opacity="', VIOLIN_STROKE_OPACITY,
 '" stroke-width="', VIOLIN_STROKE_WIDTH_PX,
 '" stroke-linejoin="round"/>'
 )

 observed <- which(!is.na(values))

 for (repeat_index in observed) {
 count <- values[repeat_index]
 # All neuronal-count dots lie on the horizontal centerline.
 point_y <- center_y

 add(
 '<circle cx="', repeat_x[repeat_index],
 '" cy="', point_y,
 '" r="', POINT_RADIUS_PX,
 '" fill="', color,
 '" fill-opacity="', POINT_OPACITY, '"/>'
 )

 if (SHOW_COUNT_LABELS) {
 add(
 '<text x="', repeat_x[repeat_index],
 '" y="', center_y - 7,
 '" text-anchor="middle" font-family="Arial, sans-serif" ',
 'font-size="8" fill="#2F353C">',
 count,
 '</text>'
 )
 }
 }
 }
 }

 add('</svg>')

 writeLines(svg, output_file, useBytes = TRUE)

 if (!file.exists(output_file)) {
 stop("The SVG was not created: ", output_file)
 }

 # Hard validation: the finished SVG must contain no animal-ID prefix.
 saved_text <- paste(readLines(output_file, warn = FALSE), collapse = "\n")
 forbidden_id_prefix <- paste0("Sam", "ple")

 if (grepl(forbidden_id_prefix, saved_text, fixed = TRUE)) {
 stop("Validation failed: an animal ID was found in the SVG.")
 }

 message(
 "Verified ", SCRIPT_VERSION,
 ": no animal IDs; animal spacing = ",
 ANIMAL_CENTER_GAP_PX,
 " SVG px.\n",
 "Editable SVG saved successfully at: ",
 normalizePath(output_file, winslash = "/", mustWork = TRUE)
 )

 invisible(output_file)
}


# =============================================================================
# 5. RUN, SAVE AND VISUALIZE
# =============================================================================

saved_svg <- create_horizontal_count_plot()

# This opens the saved vector figure for visual inspection.
if (interactive() && OPEN_AFTER_SAVE) {
 browseURL(saved_svg)
}