import pygame
import pymunk
import imageio
from pymunk.pygame_util import DrawOptions
import math
import pygame.gfxdraw

frame_rate = 60  # Frames per second
duration = 18  # Desired duration in seconds
num_frames = frame_rate * duration

# Define constants
WIDTH = 1350
HEIGHT = 500
BALL_RADIUS = 60
NUM_BALLS = 5
INITIAL_IMPULSE = 2000

# Initialize Pygame and create a window
pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Newton's Cradle - REACT")
clock = pygame.time.Clock()

# Initialize Pymunk and create a space
space = pymunk.Space()
space.gravity = (0, 4000)
space.damping = 0.5

# Create the balls and the text labels
balls = []
labels = []
font = pygame.font.Font(None, 120)  # Double the font size
text = "REACT"
ball_color = (36, 103, 151)  # Converted from hex #246797 to RGB
for i in range(NUM_BALLS):
    mass = 1
    moment = pymunk.moment_for_circle(mass, 0, BALL_RADIUS, (0, 0))
    body = pymunk.Body(mass, moment)
    body.position = (WIDTH // 2 - BALL_RADIUS * 4 + i * BALL_RADIUS * 2, HEIGHT // 2 + 100)
    
    shape = pymunk.Circle(body, BALL_RADIUS)
    shape.elasticity = 0.98
    shape.friction = 0.0
    space.add(body, shape)
    balls.append(shape)

    label = font.render(text[i], True, (255, 255, 255))
    labels.append(label)

# Apply force to the leftmost ball to simulate the inclined starting position
angle = 30  # Angle in degrees
force_magnitude = 500  # You can adjust the force magnitude as needed
force_x = force_magnitude * math.sin(math.radians(angle))
force_y = -force_magnitude * math.cos(math.radians(angle))
balls[4].body.apply_force_at_local_point((force_x, force_y), (0, 0))

# Create the strings (constraints) connecting the balls to the ceiling
strings = []
for ball in balls:
    anchor = pymunk.Body(body_type=pymunk.Body.STATIC)
    anchor.position = ball.body.position + (0, -BALL_RADIUS - 300)  # Adjusted position
    tether_length = (ball.body.position - anchor.position).length
    joint = pymunk.SlideJoint(anchor, ball.body, (0, 0), (0, 0), tether_length, tether_length)
    space.add(joint)
    strings.append((anchor, ball))

frames=[]

# Apply force to the leftmost ball to simulate the inclined starting position
angle = 30  # Angle in degrees
force_magnitude = 500  # You can adjust the force magnitude as needed
force_x = force_magnitude * math.sin(math.radians(angle))
force_y = -force_magnitude * math.cos(math.radians(angle))
balls[4].body.apply_force_at_local_point((force_x, force_y), (0, 0))

# Apply initial impulse to the first ball
balls[0].body.apply_impulse_at_local_point((INITIAL_IMPULSE, 0))

# Main loop
running = True
current_frame = 0
while running and current_frame < num_frames:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    screen.fill((255, 255, 255))

    # Draw the tethers
    for anchor, ball in strings:
        x1, y1 = anchor.position
        x2, y2 = ball.body.position
        tether_thickness = 3  # Change this value to control the thickness of the tethers
        for i in range(-tether_thickness // 2, tether_thickness // 2 + 1):
            for j in range(-tether_thickness // 2, tether_thickness // 2 + 1):
                pygame.draw.aaline(screen, ball_color, (x1 + i, y1 + j), (x2 + i, y2 + j))

    # Draw the balls and the text labels on the balls
    for ball, label in zip(balls, labels):
        ball_center = tuple(map(int, ball.body.position))
        pygame.gfxdraw.aacircle(screen, ball_center[0], ball_center[1], BALL_RADIUS, ball_color)
        pygame.gfxdraw.filled_circle(screen, ball_center[0], ball_center[1], BALL_RADIUS, ball_color)
        pygame.gfxdraw.aacircle(screen, ball_center[0], ball_center[1], BALL_RADIUS, ball_color)

        rect = label.get_rect()
        rect.center = ball.body.position
        screen.blit(label, rect)

    frame = pygame.Surface(screen.get_size())
    frame.blit(screen, (0, 0))
    frames.append(frame)

    pygame.display.flip()
    space.step(1 / 60.0)
    clock.tick(60)
    current_frame += 1

pygame.quit()

# Save the frames as a GIF
output_file = "newtons_cradle.gif"
with imageio.get_writer(output_file, mode='I', duration=1 / frame_rate) as writer:
    for frame in frames:
        img_array = pygame.surfarray.array3d(frame)
        img_array = img_array.swapaxes(0, 1)
        writer.append_data(img_array)
