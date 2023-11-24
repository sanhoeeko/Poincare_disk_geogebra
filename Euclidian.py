import math


def find_intersections_two_circles(circle1, circle2):
    # Unpack the circle parameters
    x1, y1, r1 = circle1
    x2, y2, r2 = circle2

    # Calculate the distance between the circle centers
    d = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5

    # Check if the circles are too far apart, too close together, or coincident
    if d > r1 + r2 or d < abs(r1 - r2) or d == 0 and r1 == r2:
        # No intersections
        return []

    # Calculate the angle between the line joining the circle centers and the x-axis
    a = math.atan2(y2 - y1, x2 - x1)

    # Calculate the distance from the first circle center to the line through the intersections
    h = ((r1 ** 2) - (r2 ** 2) + (d ** 2)) / (2 * d)

    # Calculate the coordinates of the point on that line
    x3 = x1 + h * math.cos(a)
    y3 = y1 + h * math.sin(a)

    # Calculate the distance from that point to the intersections
    k = (r1 ** 2 - h ** 2) ** 0.5

    # Calculate the coordinates of the intersections
    x4 = x3 + k * math.sin(a)
    y4 = y3 - k * math.cos(a)
    x5 = x3 - k * math.sin(a)
    y5 = y3 + k * math.cos(a)

    # Return the intersections as a list of tuples
    return [(x4, y4), (x5, y5)]


# Define a function to find the intersections of a circle and a line
def find_intersections_circle_line(circle, line):
    # Unpack the circle parameters
    x_center, y_center, radius = circle

    # Unpack the line parameters
    a, b, c = line

    # Check if the line is vertical (b = 0)
    if b == 0:
        # Find the x-coordinate of the line
        x = -c / a
        # Find the discriminant of the quadratic equation for y
        d = radius ** 2 - (x - x_center) ** 2
        # Check if the discriminant is negative, zero, or positive
        if d < 0:
            # No intersection
            return []
        elif d == 0:
            # One intersection
            y = y_center
            return [(x, y)]
        else:
            # Two intersections
            y1 = y_center + math.sqrt(d)
            y2 = y_center - math.sqrt(d)
            return [(x, y1), (x, y2)]
    else:
        # Find the slope and intercept of the line
        m = -a / b
        k = -c / b
        # Find the coefficients of the quadratic equation for x
        A = 1 + m ** 2
        B = -2 * (x_center - m * (k - y_center))
        C = (x_center ** 2 + (k - y_center) ** 2 - radius ** 2)
        # Find the discriminant of the quadratic equation for x
        d = B ** 2 - 4 * A * C
        # Check if the discriminant is negative, zero, or positive
        if d < 0:
            # No intersection
            return []
        elif d == 0:
            # One intersection
            x = -B / (2 * A)
            y = m * x + k
            return [(x, y)]
        else:
            # Two intersections
            x1 = (-B + math.sqrt(d)) / (2 * A)
            x2 = (-B - math.sqrt(d)) / (2 * A)
            y1 = m * x1 + k
            y2 = m * x2 + k
            return [(x1, y1), (x2, y2)]
