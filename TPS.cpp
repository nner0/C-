#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <fstream>

struct Point {
    double x, y;
};

double distance(const Point& a, const Point& b) {
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

std::vector<int> TSP(const std::vector<Point>& points) {
    int n = points.size();
    std::vector<bool> visited(n, false);
    std::vector<int> path;
    path.reserve(n);
    path.push_back(0);
    visited[0] = true;
    
    for (int i = 1; i < n; ++i) {
        int last = path.back();
        int next = -1;
        double best_dist = std::numeric_limits<double>::max();
        
        for (int j = 0; j < n; ++j) {
            if (!visited[j]) {
                double dist = distance(points[last], points[j]);
                if (dist < best_dist) {
                    best_dist = dist;
                    next = j;
                }
            }
        }
        
        path.push_back(next);
        visited[next] = true;
    }
    return path;
}

double path_length(const std::vector<Point>& points, const std::vector<int>& path, bool isCycle) {
    double length = 0.0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        length += distance(points[path[i]], points[path[i + 1]]);
    }
    if (isCycle) {
        length += distance(points[path.back()], points[path.front()]);
    }
    return length;
}

void two_opt(std::vector<int>& path, const std::vector<Point>& points) {
    bool improved = true;
    int n = path.size();
    while (improved) {
        improved = false;
        for (int i = 1; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double oldDist = distance(points[path[i - 1]], points[path[i]]) +
                                 distance(points[path[j]], points[path[(j + 1) % n]]);
                double newDist = distance(points[path[i - 1]], points[path[j]]) +
                                 distance(points[path[i]], points[path[(j + 1) % n]]);
                if (newDist < oldDist) {
                    std::reverse(path.begin() + i, path.begin() + j + 1);
                    improved = true;
                }
            }
        }
    }
}

void three_opt(std::vector<int>& path, const std::vector<Point>& points) {
    int n = path.size();
    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 0; i < n - 2; ++i) {
            for (int j = i + 1; j < n - 1; ++j) {
                for (int k = j + 1; k < n; ++k) {
                    std::vector<int> newPath = path;
                    std::reverse(newPath.begin() + i, newPath.begin() + j);
                    std::reverse(newPath.begin() + j, newPath.begin() + k);
                    double oldDist = path_length(points, path, true);
                    double newDist = path_length(points, newPath, true);
                    if (newDist < oldDist) {
                        path = newPath;
                        improved = true;
                    }
                }
            }
        }
    }
}

int main() {
    std::ifstream inputFile("../tsp_51_1");
    if (!inputFile) {
        std::cerr << "Error opening input file!\n";
        return 1;
    }
    
    int n;
    inputFile >> n;
    std::vector<Point> points(n);
    for (int i = 0; i < n; ++i) {
        inputFile >> points[i].x >> points[i].y;
    }
    inputFile.close();
    
    std::vector<int> path = TSP(points);
    bool isCycle = true;
    two_opt(path, points);
    three_opt(path, points);
    double length = path_length(points, path, isCycle);
    
    std::cout << length << " " << isCycle << "\n";
    for (int v : path) {
        std::cout << v + 1 << " ";
    }
    std::cout << "\n";
    
    return 0;
}