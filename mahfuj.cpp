#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <cmath>
#include <limits>
#include <algorithm> // Added for reverse()

using namespace std;

struct Node
{
    double latitude, longitude;
    bool operator<(const Node &other) const
    {
        return tie(latitude, longitude) < tie(other.latitude, other.longitude);
    }
    bool operator==(const Node &other) const
    {
        return latitude == other.latitude && longitude == other.longitude;
    }
};

struct Edge
{
    Node to;
    double weight;
};

map<Node, vector<Edge>> graph;

// Haversine formula to calculate distance between two lat-long points
constexpr double R = 6371.0; // Earth radius in km

double haversine(double lat1, double lon1, double lat2, double lon2)
{
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;
    lat1 = lat1 * M_PI / 180.0;
    lat2 = lat2 * M_PI / 180.0;
    double a = sin(dLat / 2) * sin(dLat / 2) + cos(lat1) * cos(lat2) * sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return R * c;
}

void addEdge(Node u, Node v, double weight)
{
    graph[u].push_back({v, weight});
    graph[v].push_back({u, weight});
}

void readTXT(const string &filename)
{
    ifstream file(filename);
    if (!file)
    {
        cerr << "Error opening file!" << endl;
        return;
    }
    string line;
    while (getline(file, line))
    {
        stringstream ss(line);
        string type;
        getline(ss, type, ',');
        vector<Node> nodes;
        double lat, lon;
        while (ss >> lon)
        {
            if (ss.peek() == ',')
                ss.ignore();
            if (!(ss >> lat))
                break;
            nodes.push_back({lat, lon});
            if (ss.peek() == ',')
                ss.ignore();
        }
        if (nodes.size() < 2)
            continue;
        for (size_t i = 0; i < nodes.size() - 1; i++)
        {
            double weight = haversine(nodes[i].latitude, nodes[i].longitude, nodes[i + 1].latitude, nodes[i + 1].longitude);
            addEdge(nodes[i], nodes[i + 1], weight);
        }
    }
}

void dijkstra(Node src, Node dest)
{
    map<Node, double> dist;
    map<Node, Node> parent;
    set<Node> visited;
    priority_queue<pair<double, Node>, vector<pair<double, Node>>, greater<pair<double, Node>>> pq;

    for (auto &pair : graph)
    {
        dist[pair.first] = numeric_limits<double>::infinity();
    }
    dist[src] = 0;
    pq.push({0, src});

    while (!pq.empty())
    {
        auto top = pq.top();
        pq.pop();
        double d = top.first;
        Node u = top.second;
        if (visited.count(u))
            continue;
        visited.insert(u);
        if (u == dest)
            break;

        for (auto &edge : graph[u])
        {
            Node v = edge.to;
            double weight = edge.weight;
            if (dist[u] + weight < dist[v])
            {
                dist[v] = dist[u] + weight;
                parent[v] = u;
                pq.push({dist[v], v});
            }
        }
    }

    if (dist[dest] == numeric_limits<double>::infinity())
    {
        cout << "No path found!" << endl;
        return;
    }

    cout << "Shortest path distance: " << dist[dest] << " km" << endl;
    vector<Node> path;
    for (Node at = dest; !(at == src); at = parent[at])
    {
        path.push_back(at);
    }
    path.push_back(src);
    reverse(path.begin(), path.end());
    for (auto &node : path)
    {
        cout << "(" << node.latitude << ", " << node.longitude << ") -> ";
    }
    cout << "END" << endl;
}

int main()
{
    string filename = "Roadmap-Dhaka.csv";
    readTXT(filename);
    Node src = {23.855136, 90.404772}; // Example source
    Node dest = {23.758633, 90.44006}; // Example destination
    dijkstra(src, dest);
    return 0;
}