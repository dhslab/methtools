import pandas as pd
import numpy as np
import scipy.stats as stats
from matplotlib.colors import to_rgb, to_hex


def fisher_exact_test(group1_meth, group1_unmeth, group2_meth, group2_unmeth):
    """
    Perform Fisher's Exact Test, handling potential NaN values.

    This function checks if any of the input values are NaN. If so, it returns
    NaN for both the odds ratio and p-value. Otherwise, it performs the test.

    Parameters:
    group1_meth (int): Methylated CpG count for group 1.
    group1_unmeth (int): Unmethylated CpG count for group 1.
    group2_meth (int): Methylated CpG count for group 2.
    group2_unmeth (int): Unmethylated CpG count for group 2.

    Returns:
    pd.Series: A pandas Series containing:
        - 'odds_ratio': Odds Ratio from the test (or np.nan).
        - 'p_value': p-value from the test (or np.nan).
    """
    # Check for any NaN or non-finite values in the inputs for a given row
    inputs = [group1_meth, group1_unmeth, group2_meth, group2_unmeth]
    if any(pd.isna(v) for v in inputs):
        return pd.Series({"odds_ratio": np.nan, "p_value": np.nan})

    # Ensure all values are integers for the contingency table
    table = [
        [int(group1_meth), int(group1_unmeth)],
        [int(group2_meth), int(group2_unmeth)],
    ]

    try:
        # Perform Fisher's exact test
        odds_ratio, p_value = stats.fisher_exact(table, alternative="two-sided")
        return pd.Series({"odds_ratio": odds_ratio, "p_value": p_value})
    except ValueError:
        # This can happen if a row/column in the table sums to zero
        return pd.Series({"odds_ratio": np.nan, "p_value": np.nan})


def get_scaled_corner_distance(x, y):
    """
    Calculates the scaled distance to the nearest corner for points in a unit square.

    This function takes x and y coordinates within a [0,1] x [0,1] unit square
    and computes the distance to the nearest corner. The distance is scaled such
    that a point exactly on a corner has a scaled distance of 0, and a point
    at the center of the square (0.5, 0.5) has a scaled distance of 1.

    Args:
        x (np.ndarray): A NumPy array of x-coordinates, ideally between 0 and 1.
        y (np.ndarray): A NumPy array of y-coordinates, ideally between 0 and 1.

    Returns:
        np.ndarray: A NumPy array containing the scaled distance (0 to 1) for
                    each input point to its nearest corner.

    Raises:
        ValueError: If the input arrays 'x' and 'y' do not have the same shape.
    """

    # --- 1. Input Validation ---
    # Ensure x and y are NumPy arrays for efficient vectorized operations
    x = np.asarray(x)
    y = np.asarray(y)

    if x.shape != y.shape:
        raise ValueError("Input arrays 'x' and 'y' must have the same shape.")

    if np.any((x < 0) | (x > 1)) or np.any((y < 0) | (y > 1)):
        import warnings

        warnings.warn(
            "Some x/y coordinates are outside the recommended [0,1] range.", UserWarning
        )

    # --- 2. Setup ---
    # Define the corners of the unit square
    corners = np.array(
        [
            [0, 0],  # Bottom-left
            [1, 0],  # Bottom-right
            [0, 1],  # Top-left
            [1, 1],  # Top-right
        ]
    )

    # The maximum possible "minimum distance" occurs at the center (0.5, 0.5).
    # The Euclidean distance from (0.5, 0.5) to any corner is sqrt(0.5^2 + 0.5^2) = sqrt(0.5).
    # This value is used to normalize the distances for scaling.
    max_min_dist = np.sqrt(0.5)

    # --- 3. Vectorized Distance Calculation ---
    # Create a combined array of input points (shape: [n, 2])
    points = np.stack((x, y), axis=-1)

    # Calculate distances from each point to all four corners.
    # This uses broadcasting to efficiently compute the differences.
    # `points[:, np.newaxis, :]` has shape [n, 1, 2]
    # `corners` has shape [4, 2]
    # The result `delta` has shape [n, 4, 2]
    delta = points[:, np.newaxis, :] - corners

    # Calculate Euclidean distances (shape: [n, 4])
    distances_to_corners = np.sqrt(np.sum(delta**2, axis=2))

    # Find the minimum distance for each point (shape: [n,])
    min_dist = np.min(distances_to_corners, axis=1)

    # --- 4. Scale the Distances ---
    # Scale the minimum distance by the maximum possible minimum distance.
    # We clip the values between 0 and 1 to handle potential floating-point
    # inaccuracies or points that fall outside the unit square.
    scaled_dist = np.clip(min_dist / max_min_dist, 0, 1)

    return scaled_dist


# # --- Example Usage ---

# # 1. Create a grid of points to visualize the distances
# # Using np.linspace to create evenly spaced points
# x_coords = np.linspace(0, 1, 100)
# y_coords = np.linspace(0, 1, 100)
# xx, yy = np.meshgrid(x_coords, y_coords)

# # Flatten the grid for processing
# x_flat = xx.ravel()
# y_flat = yy.ravel()

# # 2. Calculate the scaled distance for each point
# distances = get_scaled_corner_distance(x_flat, y_flat)

# # 3. Create a DataFrame for plotting
# df = pd.DataFrame({"x": x_flat, "y": y_flat, "distance": distances})

# # 4. Plot the results using lets-plot
# # We use the 'distance' for the color of the points, similar to the R example.
# # The `scale_color_gradient` creates a smooth grayscale transition.
# # The R code `gray(1 - distance)` results in white corners (dist=0) and a black
# # center (dist=1). We replicate this by setting `low='white'` and `high='black'`.
# # A gray background is used to ensure the white points are visible.
# p = (
#     ggplot(df, aes(x="x", y="y"))
#     + geom_point(
#         aes(color="distance"), shape=15, size=1.5
#     )  # shape=15 is a solid square like R's pch=15
#     + scale_color_gradient(low="white", high="black", name="Scaled\nDistance")
#     + coord_fixed(ratio=1)
#     + labs(title="Scaled Distance to Nearest Corner", x="X-axis", y="Y-axis")
# )

# p.show()


def assign_corner_color(
    x,
    y,
    corner_colors=("dodgerblue", "firebrick", "seagreen", "goldenrod"),
    far_color="gray",
    ramp_speed=1.0,
    na_color=None,
):
    """
    Assigns a color to points based on their proximity to the corners of a unit square.

    This function takes x/y coordinates within a [0,1]x[0,1] unit square and
    assigns a color to each point. The color is an interpolation between the
    color of the nearest corner and a specified 'far_color'. The function is
    vectorized for performance and handles missing (NaN) values gracefully.

    Args:
        x (array-like): A list, pandas Series, or NumPy array of x-coordinates.
        y (array-like): A list, pandas Series, or NumPy array of y-coordinates.
        corner_colors (tuple, optional): A tuple of four color strings.
            The colors correspond to the corners in this order:
            bottom-left (0,0), bottom-right (1,0), top-left (0,1),
            and top-right (1,1). Defaults to a standard set of four colors.
        far_color (str, optional): The color that points are scaled towards as
            they get farther from their closest corner. Defaults to 'gray'.
        ramp_speed (float, optional): Controls the rate of color transition.
            Values > 1 make corner colors more dominant (slower transition).
            Values < 1 make the transition to 'far_color' happen more quickly.
            Defaults to 1.0.
        na_color (str or None, optional): The color to assign to points where
            either x or y is NaN. Defaults to None, which means these points
            will typically not be plotted.

    Returns:
        np.ndarray: A NumPy array of hex color codes, one for each input point.

    Raises:
        ValueError: If input arrays 'x' and 'y' have different lengths or if
                    'corner_colors' does not contain exactly four colors.
    """
    # --- 1. Input Validation and Conversion ---
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.shape != y.shape:
        raise ValueError("Input arrays 'x' and 'y' must have the same shape.")
    if len(corner_colors) != 4:
        raise ValueError(
            "'corner_colors' must be a tuple or list of exactly four colors."
        )

    # --- 2. Handle NA values FIRST ---
    # Create a results array to store the final colors, initialized with na_color
    final_colors = np.full(x.shape, na_color, dtype=object)

    # Create a boolean mask for valid (non-NA) data
    valid_mask = ~np.isnan(x) & ~np.isnan(y)

    # If there are no valid points, return early
    if not np.any(valid_mask):
        return final_colors

    x_valid = x[valid_mask]
    y_valid = y[valid_mask]

    # --- 3. Setup for Color Calculation ---
    corners = np.array([[0, 0], [1, 0], [0, 1], [1, 1]])

    # Convert corner and far colors from names to RGB tuples (0-1 range)
    corner_rgb = np.array([to_rgb(c) for c in corner_colors])
    far_rgb = np.array(to_rgb(far_color))

    # The maximum possible "minimum distance" occurs at the center (0.5, 0.5)
    max_min_dist = np.sqrt(0.5)

    # --- 4. Vectorized Distance and Index Calculation ---
    points = np.stack((x_valid, y_valid), axis=1)

    # Calculate distances from each point to all four corners using broadcasting
    delta = points[:, np.newaxis, :] - corners
    distances_to_corners = np.sqrt(np.sum(delta**2, axis=2))

    # Find the index of the closest corner and the distance to it for each point
    closest_corner_indices = np.argmin(distances_to_corners, axis=1)
    min_distances = np.min(distances_to_corners, axis=1)

    # --- 5. Vectorized Color Interpolation ---
    # Get the RGB color of the closest corner for each point
    closest_corner_colors_rgb = corner_rgb[closest_corner_indices]

    # Calculate the interpolation proportion, adjusted by ramp_speed
    raw_proportion = np.clip(min_distances / max_min_dist, 0, 1)
    proportion = raw_proportion**ramp_speed

    # Reshape proportion for broadcasting with RGB colors
    proportion = proportion[:, np.newaxis]

    # Linearly interpolate between the corner color and the far color
    # Formula: C1 * (1 - p) + C2 * p
    interpolated_rgb = (
        closest_corner_colors_rgb * (1 - proportion) + far_rgb * proportion
    )

    # Convert the calculated RGB values back to hex codes
    calculated_colors_hex = [to_hex(c) for c in interpolated_rgb]

    # --- 6. Finalize and Return ---
    # Place the calculated hex colors into the correct positions in the final array
    final_colors[valid_mask] = calculated_colors_hex

    return final_colors


# # --- Example Usage ---

# # 1. Define colors
# corner_cols = ("dodgerblue", "firebrick", "seagreen", "goldenrod")
# far_away_color = "white"  # A dark gray, equivalent to R's "gray20"

# # 2. Create a grid of points, then introduce some NA values
# np.random.seed(123)  # for reproducibility
# x_coords = np.linspace(0, 1, 30)
# y_coords = np.linspace(0, 1, 30)
# xx, yy = np.meshgrid(x_coords, y_coords)
# point_data = pd.DataFrame({"x": xx.ravel(), "y": yy.ravel()})

# # Randomly set 10% of points to have NA in x or y
# na_indices = np.random.choice(
#     point_data.index, size=int(0.1 * len(point_data)), replace=False
# )
# point_data.loc[na_indices, "x"] = np.nan

# # 3. Calculate colors using the function
# point_data["color"] = assign_corner_color(
#     x=point_data["x"],
#     y=point_data["y"],
#     corner_colors=corner_cols,
#     far_color=far_away_color,
#     ramp_speed=0.5,
# )

# # 4. Plot the results using lets-plot
# # We use scale_color_identity() because the 'color' column already contains valid hex codes.
# p = (
#     ggplot(point_data.dropna(subset=["color"]), aes(x="x", y="y"))
#     + geom_point(aes(color="color"), size=2)
#     + scale_color_identity()
#     + coord_fixed(ratio=1)
#     + labs(
#         title="Python: Color Gradient with NA values handled", x="X-axis", y="Y-axis"
#     )
#     + theme_minimal()
# )

# p.show()


def assign_blended_corner_color(
    x,
    y,
    color1="dodgerblue",
    color2="firebrick",
    corner_strength=2.0,  # how fast corners dominate near their diagonals
    overlay_sharpness=2.0,  # how sharply the overlay ramps up toward center
    clip=True,
    na_color=None,
):
    """
    Color a [0,1]x[0,1] unit square with two base colors and a center overlay:

      • color1 is maximal at (0,0) and (1,1)  (main diagonal).
      • color2 is maximal at (0,1) and (1,0)  (anti-diagonal).
      • At (0.5,0.5) the color is the *maximum overlay* of color1 & color2
        (distinct third color), not a simple average.

    The overlay contribution increases smoothly toward the center, while the base
    color is a linear mix between color1/color2 based on proximity to the two
    diagonal corner pairs.

    Parameters
    ----------
    x, y : array-like of floats in [0,1]
        Coordinates (same shape). If clip=True, values are clipped to [0,1].
    color1, color2 : str or RGB tuple
        Base colors for the two diagonals.
    corner_strength : float >= 1
        Larger -> corners dominate more steeply near their diagonals.
    overlay_sharpness : float >= 1
        Larger -> overlay peaks more sharply near the center.
    clip : bool
        If True, clip x,y into [0,1]; else raise if outside range.
    na_color : str or None
        Color to assign where x or y is NaN. (None means "leave empty".)

    Returns
    -------
    np.ndarray of hex strings, same shape as x/y.
    """

    # --- input prep ---
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape.")
    if clip:
        x = np.clip(x, 0.0, 1.0)
        y = np.clip(y, 0.0, 1.0)
    else:
        if (x.min() < 0) or (x.max() > 1) or (y.min() < 0) or (y.max() > 1):
            raise ValueError("x and y must lie within [0,1] when clip=False.")

    out = np.full(x.shape, na_color, dtype=object)
    valid = ~np.isnan(x) & ~np.isnan(y)
    if not np.any(valid):
        return out

    xv = x[valid]
    yv = y[valid]

    # --- colors as RGB in [0,1] ---
    c1 = np.array(to_rgb(color1), dtype=float)
    c2 = np.array(to_rgb(color2), dtype=float)

    # --- distances to the TWO corner-pairs (not individual corners) ---
    # Pair 1: (0,0) and (1,1) -> color1
    d10 = np.hypot(xv - 0.0, yv - 0.0)
    d11 = np.hypot(xv - 1.0, yv - 1.0)
    d1 = np.minimum(d10, d11)  # nearest distance to pair 1; in [0,1] within unit square

    # Pair 2: (0,1) and (1,0) -> color2
    d20 = np.hypot(xv - 0.0, yv - 1.0)
    d21 = np.hypot(xv - 1.0, yv - 0.0)
    d2 = np.minimum(d20, d21)  # nearest distance to pair 2; in [0,1] within unit square

    # Convert distances to *corner* weights (bigger near their diagonals)
    # Normalize by 1.0: at the opposite pair's corners, min-distance is exactly 1.
    w1_raw = np.power(1.0 - np.clip(d1, 0.0, 1.0), corner_strength)
    w2_raw = np.power(1.0 - np.clip(d2, 0.0, 1.0), corner_strength)

    s = w1_raw + w2_raw
    s = np.where(s == 0, 1.0, s)  # safety
    w1 = w1_raw / s
    w2 = w2_raw / s  # ensures w1 + w2 = 1

    # --- overlay amount: 0 at corners, 1 at center, smooth in between ---
    # Using 4*w1*w2 ∈ [0,1], max=1 when w1=w2=0.5 (center line), ~0 at corners.
    t = np.power(4.0 * w1 * w2, overlay_sharpness)

    # --- base linear mix (without overlay) ---
    base = (w1[:, None] * c1) + (w2[:, None] * c2)

    # --- symmetric OVERLAY blend of the two *fixed* colors ---
    # Standard overlay per channel; we average overlay(c1,c2) and overlay(c2,c1)
    # to make it commutative.
    def _overlay(a, b):
        # a, b are in [0,1]; vectorized over channels
        return np.where(a < 0.5, 2.0 * a * b, 1.0 - 2.0 * (1.0 - a) * (1.0 - b))

    O1 = _overlay(c1, c2)
    O2 = _overlay(c2, c1)
    overlay_color = 0.5 * (O1 + O2)  # shape (3,)

    # --- combine: lerp between base mix and overlay color by t ---
    rgb = ((1.0 - t)[:, None] * base) + (t[:, None] * overlay_color)

    # to hex
    colors_hex = np.array([to_hex(c) for c in rgb], dtype=object)
    out[valid] = colors_hex
    return out
