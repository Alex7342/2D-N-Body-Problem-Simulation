#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <iostream>
#include <thread>


// Constants
const double G = 6.67430e-11;   // Gravitational constant
const double TIME_STEP = 1;   // Time step in seconds
const int STEPS = 4000;        // Number of steps to simulate
const int WINDOW_WIDTH = 1600;  // The width of the window
const int WINDOW_HEIGHT = 900;  // The height of the window
const double SCALE = 0.1;       // Scale factor for visualization
const int NUMBER_OF_THREADS = 8;  // Number of threads

// Struct to represent a body in the simulation
struct Body {
    double mass;
    double x, y;        // Position
    double vx, vy;      // Velocity
    sf::CircleShape shape; // Visual representation
    sf::Color color;    // Color of the object
};

// Function to generate a random body structure
Body generateRandomBody(int minX, int maxX, int minY, int maxY) {
    double mass = (1 + rand() % 9) * pow(10, 7 + rand() % 5);
    double x = minX + rand() % (maxX - minX);
    double y = minY + rand() % (maxY - minY);
    double xv = 1.0 * (rand() % 4) / 10.0;
    double yv = 1.0 * (rand() % 4) / 10.0;
    return Body{ mass, x, y, xv, yv, sf::CircleShape(3), sf::Color::White };
}

Body generateRandomStar(int minX, int maxX, int minY, int maxY) {
    double mass = (1 + rand() % 9) * pow(10, 13);
    double x = minX + rand() % (maxX - minX);
    double y = minY + rand() % (maxY - minY);
    double xv = 1.0 * (rand() % 4) / 10.0;
    double yv = 1.0 * (rand() % 4) / 10.0;
    return Body{ mass, x, y, xv, yv, sf::CircleShape(6), sf::Color::Yellow };
}

// Function to compute the distance between two bodies
double computeDistance(const Body& a, const Body& b) {
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    return std::sqrt(dx * dx + dy * dy);
}

// Function to compute the gravitational force between two bodies
void computeForce(const Body& a, const Body& b, double& fx, double& fy) {
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double distance = computeDistance(a, b);
    double force = (G * a.mass * b.mass) / (distance * distance + 1e-10);
    fx = force * (dx / distance);
    fy = force * (dy / distance);
}

void threadFunction(std::vector<Body>& bodies, const int begin, const int end) {
    for (size_t i = begin; i < end; ++i) {
        Body& body = bodies[i];

        // Calculate gravitational force
        std::pair<double, double> force{ 0, 0 };
        for (size_t j = 0; j < bodies.size(); ++j) {
            if (i != j) {
                Body& otherBody = bodies[j];
                double fx, fy;
                computeForce(body, otherBody, fx, fy);
                force.first += fx;
                force.second += fy;
            }
        }

        // Update positions and velocities
        double ax = force.first / body.mass;
        double ay = force.second / body.mass;

        body.vx += ax * TIME_STEP;
        body.vy += ay * TIME_STEP;

        body.x += body.vx * TIME_STEP;
        body.y += body.vy * TIME_STEP;

        // Update the visual representation
        body.shape.setPosition(
            WINDOW_WIDTH / 2 + body.x * SCALE,
            WINDOW_HEIGHT / 2 + body.y * SCALE
        );

        body.shape.setFillColor(body.color);
    }
}

// Main simulation and visualization function
void simulateAndVisualize(std::vector<Body>& bodies) {
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "N-Body Simulation");
    int currentStep = 0;

    std::chrono::system_clock::time_point startTime = std::chrono::system_clock::now();

    while (window.isOpen() && currentStep < STEPS) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        std::vector<std::thread> threads;
        for (int i = 0; i < NUMBER_OF_THREADS; i++) {
            int begin = (i * bodies.size()) / NUMBER_OF_THREADS;
            int end = ((i + 1) * bodies.size()) / NUMBER_OF_THREADS;
            if (i == NUMBER_OF_THREADS - 1) end = bodies.size();
            threads.push_back(std::thread(std::ref(threadFunction), std::ref(bodies), begin, end));
        }

        for (int i = 0; i < NUMBER_OF_THREADS; i++)
            threads[i].join();

        // Render
        window.clear(sf::Color::Black);
        for (const auto& body : bodies) {
            window.draw(body.shape);
        }
        window.display();

        // Output the state of the system
        /*for (const auto& body : bodies) {
            std::cout << "Body at (" << body.x << ", " << body.y << ")\n";
        }
        std::cout << "------------------------\n";*/

        currentStep++;
    }

    std::chrono::system_clock::time_point stopTime = std::chrono::system_clock::now();
    std::cout << "Time to simulate " << currentStep << " steps of " << TIME_STEP << " seconds each for "
        << bodies.size() << " bodies, using " << NUMBER_OF_THREADS << " threads: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count() << "ms\n";
}

int main() {
    srand(time(0));

    // Initialize bodies
    std::vector<Body> bodies;
    for (int i = 0; i < 500; i++)
        bodies.push_back(generateRandomBody(-7000, 7000, -4000, 4000));
    for (int i = 0; i < 3; i++)
        bodies.push_back(generateRandomStar(-7000, 7000, -4000, 4000));

    // Massive central body
    //bodies.push_back({ 5.0e13, 0, 0, 0, 0, sf::CircleShape(14), sf::Color::Yellow });

    // Set colors and origins for visualization
    for (auto& body : bodies) {
        body.shape.setFillColor(sf::Color::White);
        body.shape.setOrigin(body.shape.getRadius(), body.shape.getRadius());
    }

    // Run simulation with visualization    
    simulateAndVisualize(bodies);

    return 0;
}
