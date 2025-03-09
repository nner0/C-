#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>

struct Node {
    int level, value, weight;
    double bound;
};

struct Item {
    int value, weight;
    double ratio;
};

bool cmp(const Item &a, const Item &b) {
    return a.ratio > b.ratio;
}

double bound(const Node &u, int n, int W, const std::vector<Item> &items) {
    if (u.weight >= W) return 0;
    
    double profit_bound = u.value;
    int j = u.level + 1;
    int total_weight = u.weight;
    
    while (j < n && total_weight + items[j].weight <= W) {
        total_weight += items[j].weight;
        profit_bound += items[j].value;
        j++;
    }
    
    if (j < n)
        profit_bound += (double)(W - total_weight) * items[j].ratio;
    
    return profit_bound;
}

struct CompareBound {
    bool operator()(const Node &a, const Node &b) {
        return a.bound < b.bound; //  max-heap по bound
    }
};

std::pair<int, int> knapsack(int W, std::vector<Item> &items, int n) {
    std::sort(items.begin(), items.end(), cmp);
    
    std::priority_queue<Node, std::vector<Node>, CompareBound> PQ;
    Node u, v;
    
    u.level = -1;
    u.value = u.weight = 0;
    u.bound = bound(u, n, W, items);
    PQ.push(u);
    
    int max_profit = 0;
    int total_weight = 0;
    
    while (!PQ.empty()) {
        u = PQ.top();
        PQ.pop();
        
        if (u.level + 1 >= n) continue; // чтобы не выйти за границы        
        v.level = u.level + 1;
        
        v.weight = u.weight + items[v.level].weight;
        v.value = u.value + items[v.level].value;
        
        if (v.weight <= W && v.value > max_profit) {
            max_profit = v.value;
            total_weight = v.weight;
        }
        
        v.bound = bound(v, n, W, items);
        if (v.bound > max_profit)
            PQ.push(v);

        v.weight = u.weight;
        v.value = u.value;
        v.bound = bound(v, n, W, items);
        
        if (v.bound > max_profit)
            PQ.push(v);
    }
    
    return {max_profit, total_weight};
}

int main() {
    std::ifstream input("../ks_19_0");
    if (!input) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }
    
    int N, W;
    input >> N >> W;
    std::vector<Item> items(N);
    
    for (int i = 0; i < N; i++) {
        input >> items[i].value >> items[i].weight;
        items[i].ratio = (double)items[i].value / items[i].weight;
    }
    
    input.close();
    
    if (N == 0 || W == 0) {
        std::cout << "Maximum profit: 0" << std::endl;
        std::cout << "Total weight: 0" << std::endl;
        return 0;
    }
    
    auto result = knapsack(W, items, N);
    
    std::cout << "Maximum profit: " << result.first << std::endl;
    std::cout << "Total weight: " << result.second << std::endl;
    
    return 0;
}
