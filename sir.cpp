#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <queue>
#include <limits>
#include <algorithm>
#include <iomanip>

using namespace std;

// Structure to represent a node (location)
struct Node
{
    double latitude;
    double longitude;
};

// Structure to represent an edge between nodes
struct Edge
{
    int source;
    int destination;
    double distance;
    double altitude;
    string type;            // "road" or "metro" for the edge
    string sourceName;      // For Metro, the start station
    string destinationName; // For Metro, the end station
};

struct DelimitedString
{
    string value;
    char delimiter; // delimiter that ends the value

    DelimitedString(char delim = ',') : delimiter(delim)
    {
    }

    friend istream &operator>>(istream &is, DelimitedString &output)
    {
        getline(is, output.value, output.delimiter);
        return is;
    }
};

// Function to calculate the distance between two coordinates (Haversine formula)
double calculateDistance(double lat1, double lon1, double lat2, double lon2)
{
    double R = 6371; // Radius of the Earth in kilometers
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;
    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(lat1 * M_PI / 180.0) * cos(lat2 * M_PI / 180.0) *
                   sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return R * c;
}

// Function to parse coordinates from a stringstream
vector<pair<double, double>> parseCoordinates(stringstream &ss)
{
    vector<pair<double, double>> coordinates;
    double lat, lon;
    char comma;
    DelimitedString latStr(','), lonStr(',');

    while (ss.good())
    {
        ss >> latStr >> lonStr;

        if (!ss.fail())
        {
            try
            {
                lat = stod(latStr.value);
                lon = stod(lonStr.value);
                coordinates.push_back({lat, lon});
            }
            catch (...)
            {
                ss.clear(); // Clear fail state after stod() failure
                break;
            }
        }
        else
        {
            break;
        }
    }
    return coordinates;
}
// Function to read DhakaStreet data from sample.csv
void readDhakaStreet(const string &filename, vector<Node> &nodes, vector<Edge> &edges,
                     map<pair<double, double>, int> &nodeIndex, map<int, vector<int>> &adjacencyList)
{
    ifstream file(filename);
    string line;

    if (file.is_open())
    {
        while (getline(file, line))
        {
            if (line.empty() || line[0] == '#')
            {
                cout << "Empty line detected" << endl;
                continue;
            }
            stringstream ss(line);

            DelimitedString dataType(','), nextValue(','), altitudeStr(','), distanceStr(',');
            vector<pair<double, double>> coordinates;
            double altitude = 0;
            double distance = 0;

            // Parse the line
            ss >> dataType;

            // Read coordinate pairs and altitude/distance
            while (ss.good())
            {
                // Attempt to read latitude and longitude
                DelimitedString latStr(','), lonStr(',');

                ss >> latStr; // Read latitude

                if (ss.fail())
                {
                    // No more coordinate pairs; attempt to read altitude and distance
                    try
                    {
                        altitude = stod(latStr.value);
                        ss >> distanceStr;
                        distance = stod(distanceStr.value);
                    }
                    catch (...)
                    {
                        // If the string doesn't have altitude and distance, break
                        break;
                    }
                    break;
                }

                ss >> lonStr; // Read longitude
                if (ss.fail())
                    break;

                try
                {
                    double lat = stod(latStr.value);
                    double lonValue = stod(lonStr.value);
                    coordinates.push_back({lat, lonValue});
                }
                catch (...)
                {
                    // If the string doesn't have lat long, break
                    break;
                }
            }
            // Function to add a node, only add if it doesn't exist.
            auto addNode = [&](double latitude, double longitude)
            {
                pair<double, double> coords = {latitude, longitude};
                if (nodeIndex.find(coords) == nodeIndex.end())
                {
                    Node newNode = {latitude, longitude};
                    nodes.push_back(newNode);
                    int index = nodes.size() - 1;
                    nodeIndex[coords] = index;
                    return index;
                }
                else
                {
                    return nodeIndex[coords];
                }
            };

            // Create edges between consecutive coordinate pairs
            for (size_t i = 0; i < coordinates.size() - 1; ++i)
            {
                int sourceNodeIndex = addNode(coordinates[i].first, coordinates[i].second);

                int destNodeIndex = addNode(coordinates[i + 1].first, coordinates[i + 1].second);

                Edge edge = {
                    sourceNodeIndex, destNodeIndex,
                    distance,
                    altitude, "road",
                    "",
                    ""};

                edges.push_back(edge);
                adjacencyList[sourceNodeIndex].push_back(destNodeIndex);
            }
        }
        file.close();
    }
    else
    {
        cerr << "Unable to open file " << filename << endl;
    }
}

// Function to read DhakaMetroRail data from Routemap-DhakaMetroRail.csv
void readDhakaMetroRail(const string &filename, vector<Node> &nodes, vector<Edge> &edges, map<pair<double, double>, int> &nodeIndex, map<int, vector<int>> &adjacencyList)
{
    ifstream metrofile(filename);
    string metroline;
    string transportType, startName, endName;

    if (metrofile.is_open())
    {
        while (getline(metrofile, metroline))
        {
            // Skip blank lines and lines starting with comments
            if (metroline.empty() || metroline[0] == '#')
                continue;

            stringstream ss(metroline);
            DelimitedString transportTypeStr(','), latStr(','), lonStr(','), startNameStr(','), endNameStr(',');
            vector<pair<double, double>> coordinates;

            // Parse the line
            ss >> transportTypeStr >> latStr;

            // Read coordinate pairs until no more can be read
            while (ss.good())
            {
                DelimitedString lon(',');
                ss >> lon; // Try to read the Longitude

                if (!ss.fail()) // if it was read successfully.
                {
                    double lat = 0;
                    double lonValue = 0;
                    try
                    {
                        lat = stod(latStr.value);
                        lonValue = stod(lon.value);
                    }
                    catch (...)
                    { // catch all the exceptions
                        // break from the loop to read the next csv line
                        break;
                    }

                    coordinates.push_back({lat, lonValue});

                    ss >> latStr; // load next Latitude or the startName

                    if (ss.fail()) // startName failed. we must break;
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }

            try
            {
                ss >> startNameStr >> endNameStr;
            }
            catch (...)
            {
                // cout <<"Exception occured reading start and end name"<< endl;
            }

            // Function to add a node, only add if it doesn't exist.
            auto addNode = [&](double latitude, double longitude)
            {
                pair<double, double> coords = {latitude, longitude};
                if (nodeIndex.find(coords) == nodeIndex.end())
                {
                    Node newNode = {latitude, longitude};
                    nodes.push_back(newNode);
                    int index = nodes.size() - 1;
                    nodeIndex[coords] = index;
                    return index;
                }
                else
                {
                    return nodeIndex[coords];
                }
            };
            double altitude = 0;

            if (coordinates.size() >= 2)
            {
                int sourceNodeIndex = addNode(coordinates[0].first, coordinates[0].second);
                int destNodeIndex = addNode(coordinates.back().first, coordinates.back().second);
                double total_distance = 0;

                for (size_t i = 0; i < coordinates.size() - 1; ++i)
                {
                    total_distance += calculateDistance(coordinates[i].first, coordinates[i].second, coordinates[i + 1].first, coordinates[i + 1].second);
                }

                Edge edge = {sourceNodeIndex, destNodeIndex, total_distance, altitude, "metro", startNameStr.value, endNameStr.value};
                edges.push_back(edge);
                if (adjacencyList.find(sourceNodeIndex) == adjacencyList.end())
                {
                    adjacencyList[sourceNodeIndex] = vector<int>(); // Initialize the vector if it doesn't exist
                }
                adjacencyList[sourceNodeIndex].push_back(destNodeIndex);
            }
        }
        metrofile.close();
    }
    else
    {
        cerr << "Unable to open file " << filename << endl;
    }
}

// Dijkstra's Algorithm function
pair<vector<int>, double> dijkstra(int sourceNode, int destNode, const vector<Node> &nodes, const vector<Edge> &edges,
                                   const map<int, vector<int>> &adjacencyList)
{
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq; // distance , node
    vector<double> dist(nodes.size(), numeric_limits<double>::infinity());
    vector<int> prev(nodes.size(), -1); // Store previous node in the shortest path

    dist[sourceNode] = 0;
    pq.push({0, sourceNode});

    while (!pq.empty())
    {
        double d = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        if (d > dist[u])
            continue; // Skip if we have already processed a shorter path to u

        // Iterate through all adjacent nodes of u
        if (adjacencyList.find(u) != adjacencyList.end())
        {
            for (int v : adjacencyList.at(u))
            {
                // Find the edge between u and v;
                double weight = numeric_limits<double>::infinity();
                for (const auto &edge : edges)
                {
                    if (edge.source == u && edge.destination == v)
                    {
                        weight = edge.distance;
                        break;
                    }
                }

                if (dist[v] > dist[u] + weight)
                {
                    dist[v] = dist[u] + weight;
                    prev[v] = u;
                    pq.push({dist[v], v});
                }
            }
        }
    }

    // Reconstruct the path
    vector<int> path;
    if (dist[destNode] != numeric_limits<double>::infinity())
    {
        int current = destNode;
        while (current != -1)
        {
            path.push_back(current);
            current = prev[current];
        }
        reverse(path.begin(), path.end());
    }

    return {path, dist[destNode]};
}

int main(int argc, char *argv[])
{
    std::cout << std::fixed << std::setprecision(6);

    // Data structures to store the graph
    vector<Node> nodes;
    vector<Edge> edges;
    map<pair<double, double>, int> nodeIndex; // Map coordinates to node index
    map<int, vector<int>> adjacencyList;      // adjacency list to keep track of which node goes to which nodes

    // Function to add a node, only add if it doesn't exist.
    auto addNode = [&](double latitude, double longitude)
    {
        pair<double, double> coords = {latitude, longitude};
        if (nodeIndex.find(coords) == nodeIndex.end())
        {
            Node newNode = {latitude, longitude};
            nodes.push_back(newNode);
            int index = nodes.size() - 1;
            nodeIndex[coords] = index;
            return index;
        }
        else
        {
            return nodeIndex[coords];
        }
    };

    // Read DhakaStreet and DhakaMetroRail data
    readDhakaStreet(argv[1], nodes, edges, nodeIndex,
                    adjacencyList);
    readDhakaMetroRail(argv[2], nodes, edges, nodeIndex,
                       adjacencyList);

    // Print the graph information (for verification)
    cout << "Number of nodes: " << nodes.size() << endl;
    cout << "Number of edges: " << edges.size() << endl;

    cout << "\nEdges:" << endl;
    for (const auto &edge : edges)
    {
        cout << "Type: " << edge.type << ", Source: (" << nodes[edge.source].latitude << ", "
             << nodes[edge.source].longitude << ")";

        if (edge.type == "metro")
        {
            cout << " [" << edge.sourceName << "]"; // Only print station names for Metro
            cout << ", Destination: (" << nodes[edge.destination].latitude << ", " << nodes[edge.destination].longitude
                 << ") [" << edge.destinationName << "]";
        }
        else
        {
            cout << ", Destination: (" << nodes[edge.destination].latitude << ", " << nodes[edge.destination].longitude
                 << ")";
        }

        cout << ", Distance: " << edge.distance << " km, Altitude: " << edge.altitude << endl;
    }
    // Now get Source and destination node.
    double sourceLat, sourceLon, destLat, destLon;

    sourceLat = 90.439973;
    sourceLon = 23.741813;

    destLat = 90.439341;
    destLon = 23.742218;

    // cout << "Enter the Source Latitude" << endl;
    //
    // // cin >> sourceLat;
    // cout << "Enter the Source Longitude" << endl;
    // // cin >> sourceLon;
    //
    // cout << "Enter the Destination Latitude" << endl;
    // // cin >> destLat;
    //
    // cout << "Enter the Destination Longitude" << endl;
    // // cin >> destLon;

    int sourceNode = addNode(sourceLat, sourceLon);
    int destNode = addNode(destLat, destLon);

    cout << "Source Node: " << sourceNode << endl;
    cout << "Dest Node: " << destNode << endl;

    // Run Dijkstra's algorithm
    auto [path, totalDistance] = dijkstra(sourceNode, destNode, nodes, edges, adjacencyList);

    // Print the shortest path from source to destination
    cout << "\nShortest Path from (" << nodes[sourceNode].latitude << ", " << nodes[sourceNode].longitude << ") to ("
         << nodes[destNode].latitude << ", " << nodes[destNode].longitude << "):" << endl;

    if (path.empty())
    {
        cout << "No path exists." << endl;
    }
    else
    {
        cout << "Path: ";
        for (size_t i = 0; i < path.size(); ++i)
        {
            cout << "(" << nodes[path[i]].latitude << ", " << nodes[path[i]].longitude << ")";
            if (i < path.size() - 1)
            {
                // Find Edge type
                string edgeType = "None";
                string sourceName = "";
                string destinationName = "";

                for (const auto &edge : edges)
                {
                    if (edge.source == path[i] && edge.destination == path[i + 1])
                    {
                        edgeType = edge.type;
                        sourceName = edge.sourceName;
                        destinationName = edge.destinationName;
                        break;
                    }
                }
                cout << " --(" << edgeType << ")--> ";
                if (edgeType == "metro")
                {
                    cout << "- [" << sourceName << " to " << destinationName << "]";
                }
                cout << ")--> ";
            }
        }
        cout << endl;

        cout << "Total Distance: " << totalDistance << " km" << endl;
    }

    return 0;
}
