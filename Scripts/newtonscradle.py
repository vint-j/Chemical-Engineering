import pygame
import pymunk
import imageio
from pymunk.pygame_util import DrawOptions

# Define constants
WIDTH = 600
HEIGHT = 200
BALL_RADIUS = 30
NUM_BALLS = 5
INITIAL_IMPULSE = 500

# Initialize Pygame and create a window
pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Newton's Cradle - REACT")
clock = pygame.time.Clock()

# Initialize Pymunk and create a space
space = pymunk.Space()
space.gravity = (0, 3000)

# Create the balls and the text labels
balls = []
labels = []
font = pygame.font.Font(None, 60)
text = "REACT"
ball_color = (36, 103, 151)  # Converted from hex #246797 to RGB
for i in range(NUM_BALLS):
    mass = 1
    moment = pymunk.moment_for_circle(mass, 0, BALL_RADIUS, (0, 0))
    body = pymunk.Body(mass, moment)
    body.position = (WIDTH // 2 - BALL_RADIUS * 4 + i * BALL_RADIUS * 2, HEIGHT // 2 + 50)  # Adjusted position
    shape = pymunk.Circle(body, BALL_RADIUS)
    shape.elasticity = 1.0
    shape.friction = 0.0
    space.add(body, shape)
    balls.append(shape)

    label = font.render(text[i], True, (255, 255, 255))
    labels.append(label)
    
# Create the strings (constraints) connecting the balls to the ceiling
strings = []
for ball in balls:
    anchor = pymunk.Body(body_type=pymunk.Body.STATIC)
    anchor.position = ball.body.position + (0, -BALL_RADIUS - 150)  # Adjusted position
    tether_length = (ball.body.position - anchor.position).length
    joint = pymunk.SlideJoint(anchor, ball.body, (0, 0), (0, 0), tether_length, tether_length)
    space.add(joint)
    strings.append((anchor, ball))


# Draw the balls and the text labels on the balls
for ball, label in zip(balls, labels):
    ball_center = tuple(map(int, ball.body.position))
    pygame.draw.circle(screen, ball_color, ball_center, BALL_RADIUS)  # Use ball_color instead of (255, 255, 255)


# Apply initial impulse to the first ball
balls[0].body.apply_impulse_at_local_point((INITIAL_IMPULSE, 0))

# Main loop
running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    screen.fill((255, 255, 255))

    # Draw the tethers
    for anchor, ball in strings:
        pygame.draw.line(screen, (36, 103, 151), anchor.position, ball.body.position, 2)

    # Draw the balls and the text labels on the balls
    for ball, label in zip(balls, labels):
        ball_center = tuple(map(int, ball.body.position))
        pygame.draw.circle(screen, (36, 103, 151), ball_center, BALL_RADIUS)

        rect = label.get_rect()
        rect.center = ball.body.position
        screen.blit(label, rect)

    pygame.display.flip()
    space.step(1 / 60.0)
    clock.tick(60)


pygame.quit()
