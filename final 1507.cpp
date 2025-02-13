#include <bits/stdc++.h>
using namespace std;

#define INT unsigned int
#define INF (INT)1e8

// check please: https://www.google.com/maps/d/u/0/viewer?mid=1dOfrWHkEHyWvLkNnuNhVij5PsH6PPZQ&ll=23.777787056022266%2C90.39111249999999&z=13

typedef struct Node
{
    double lat;
    double lng;
} Node;

typedef struct Edge
{
    int source;      // eta ami vabtam NODE e rakha lagbe
    int destination; // same as above
    double distance; // age etar name ami weight rakhsilam
    double altitude;

    string type;            // "road" or "metro" for the edge
    string sourceName;      // For Metro, the start station
    string destinationName; // For Metro, the end station
} Edge;

/**
 * Struct to represent a delimited string
 * This is used to read a string from a stream until a delimiter is encountered. Here delimiter is a comma.
 */
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

    // sir, ekhane tahole comma separate kore just value ta read kore rakha hobe. eta kore amra ekta string read korte pari.
};

// Function to calculate the distance between two coordinates (Haversine formula)
double calculateDistance(double lat1, double lon1, double lat2, double lon2)
{
    double R = 6371; // Radius of the Earth in kilometers
    // basically higher math e amra s = r*theta * pi / 180 use kore Dhaka to Rajsahi distance type math krsilam. erokom kichu ami bujtesi.
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;
    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(lat1 * M_PI / 180.0) * cos(lat2 * M_PI / 180.0) *
                   sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a)); // eta bujtesi na.
    return R * c;
}

// dijkstra's algorithm for this particular problem
/**
 * @param srcNode Source node index. Basically ei kaj ta amar mathay astesilo na je ami double type er long-lat ke kivabe int type er ekta node hisebe store korbo
 * @param destNode Destination node index. Same as above.
 * @param nodes Set of nodes. eta ami ageo banaisilam
 * @param edges Set of edges. eta ami ageo banaisilam. But format e ektu difference silo.
 * @param adjList Adjacency list representation of the graph. eta ami ageo banaisilam. but ami SET use korsilam. pore ar hisab milate pari ni.
 */
pair<vector<int>, double> dijkstra(int srcNode, int destNode, const vector<Node> &nodes, const vector<Edge> &edges, const map<int, vector<int>> &adjList)
{
    // Nicher code ta normal ekta dijkstra algorithm er graph
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq; // distance , node

    vector<double> dist(nodes.size(), numeric_limits<double>::infinity());
    vector<int> prev(nodes.size(), -1); // Store previous node in the shortest path. Also ekhane sob node ke -1 kore init kora hoise.

    dist[srcNode] = 0;
    pq.push({0, srcNode});

    while (!pq.empty())
    {
        // ekhane 'd' hocche node er distance, 'u' hocche node er index
        double d = pq.top().first;
        int u = pq.top().second;

        pq.pop();

        if (d > dist[u]) // distance already kom, so bypass krbo
            continue;
        else
        {
            // check korbo node gular distance
            if (adjList.find(u) != adjList.end())
            {
                for (int v : adjList.at(u))
                {
                    // Find the edge between u and v;
                    double cost = numeric_limits<double>::infinity(); // first e cost er value INF kore nilam
                    for (const auto &edge : edges)
                    {
                        if (edge.source == u && edge.destination == v) // basically srcnode r destnode er moddhe edge ache kina check korlam
                        {
                            cost = edge.distance;
                            break;
                        }
                    }

                    // nicher part ta normal dijkstra er code
                    if (dist[v] > dist[u] + cost)
                    {
                        dist[v] = dist[u] + cost;
                        prev[v] = u;
                        pq.push({dist[v], v});
                    }
                }
            }
        }
    }

    // ekhon path ta print/store korbo
    vector<int> path;
    if (dist[destNode] != numeric_limits<double>::infinity())
    {
        INT curNode = destNode;
        while (curNode != -1)
        {
            path.push_back(curNode);
            curNode = prev[curNode];
        }
        reverse(path.begin(), path.end());
    }
    return {path, dist[destNode]};

    // sir, niche amar same code ta ache. just ami store na kore print kortam normally int type node er case e.

    // rest of the code was my OLD CODE for Dijkstra's algorithm
    // vector<INT> dist(n, INF);
    // vector<INT> parent(n, -1);
    // vector<bool> visited(n, false);

    // dist[0] = 0;
    // for (INT i = 0; i < n; i++)
    // {
    //     INT min_dist = INF;
    //     INT min_index = -1;
    //     for (INT j = 0; j < n; j++)
    //     {
    //         if (!visited[j] && dist[j] < min_dist)
    //         {
    //             min_dist = dist[j];
    //             min_index = j;
    //         }
    //     }

    //     visited[min_index] = true;

    //     for (INT j = 0; j < n; j++)
    //     {
    //         if (!visited[j] && graph[min_index][j] > 0 && dist[j] > dist[min_index] + graph[min_index][j])
    //         {
    //             dist[j] = dist[min_index] + graph[min_index][j];
    //             parent[j] = min_index;
    //         }
    //     }
    // }

    // // Print the shortest path
    // for (INT i = 0; i < n; i++)
    // {
    //     cout << "Shortest path from 1 to " << i + 1 << " is: ";
    //     INT j = i;
    //     while (j != -1)
    //     {
    //         cout << j + 1 << " ";
    //         j = parent[j];
    //     }
    //     cout << "with cost: " << dist[i] << endl;
    // }
}

// Function to read DhakaMetroRail data from Routemap-DhakaMetroRail.csv
void readDhakaRoute(const string &filename, string edgeType, vector<Node> &nodes, vector<Edge> &edges, map<pair<double, double>, int> &nodeIndex, map<int, vector<int>> &adjacencyList)
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

                Edge edge = {sourceNodeIndex, destNodeIndex, total_distance, altitude, edgeType, startNameStr.value, endNameStr.value};
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

// Function to read DhakaStreet data from Roadmap-Dhaka.csv
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
                // coordinates[i] and coordinates[i + 1] are the source and destination coordinates

                // coordinates[i].first, coordinates[i].second are the lat and long of the source node
                int srcNodeId = addNode(coordinates[i].first, coordinates[i].second);

                int destNodeId = addNode(coordinates[i + 1].first, coordinates[i + 1].second);

                Edge edge = {
                    srcNodeId, destNodeId,
                    distance,
                    altitude, "road",
                    "",
                    ""};

                edges.push_back(edge);
                adjacencyList[srcNodeId].push_back(destNodeId);
            }
        }
        file.close();
    }
    else
    {
        cerr << "Unable to open file " << filename << endl;
    }
}
// Function to add a node, only add if it doesn't exist.
auto addNode(double latitude, double longitude, map<pair<double, double>, int> &nodeId, vector<Node> &nodes)
{
    pair<double, double> coords = {latitude, longitude};
    if (nodeId.find(coords) == nodeId.end())
    {
        Node newNode = {latitude, longitude};
        nodes.push_back(newNode);
        int index = nodes.size() - 1;
        nodeId[coords] = index;
        return index;
    }
    else
    {
        return nodeId[coords];
    }
};

void createKMLFile(const vector<Node> &nodes, const vector<Edge> &edges, const vector<int> &path)
{
    ofstream kmlFile("route.kml");
    if (kmlFile.is_open())
    {
        kmlFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        kmlFile << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
        kmlFile << "<Document>\n";
        kmlFile << "<name>Route</name>\n";
        kmlFile << "<Placemark>\n";
        kmlFile << "<name>bsse 1507</name>\n";
        kmlFile << "<Style>\n";
        kmlFile << "<LineStyle>\n";
        kmlFile << "<color>ffff00ff</color>\n"; // Blue color in ABGR format
        kmlFile << "<width>4</width>\n";
        kmlFile << "</LineStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<LineString>\n";
        kmlFile << "<tessellate>1</tessellate>\n";
        kmlFile << "<coordinates>\n";

        for (int i = 0; i < path.size(); ++i)
        {
            kmlFile << fixed << setprecision(6) << nodes[path[i]].lat << "," << nodes[path[i]].lng;
            kmlFile << fixed << setprecision(0);
            bool isAltitude = false;
            // check from edge and put altitude and distance
            if (i < path.size() - 1)
            {
                for (const auto &edge : edges)
                {
                    if (edge.source == path[i] && edge.destination == path[i + 1])
                    {
                        kmlFile << "," << edge.altitude << endl;
                        isAltitude = true;
                        break;
                    }
                }
            }
            if (!isAltitude)
            {
                kmlFile << ",0" << endl;
            }
        }

        kmlFile << "</coordinates>\n";
        kmlFile << "</LineString>\n";
        kmlFile << "</Placemark>\n";
        kmlFile << "</Document>\n";
        kmlFile << "</kml>\n";
        kmlFile.close();
    }
    else
    {
        cerr << "Unable to open file output.kml" << endl;
    }
}

// main function
int main(const int argc, const char *argv[])
{
    // set precision upto 6 digit
    cout << fixed << setprecision(6);
    // arguement theke file nicchi, tai argc == 2 hote hobe
    if (argc != 5)
    {
        cerr << "Usage 5 arg: " << argv[0] << " <data-set-file>" << "<route-data-set-files>(total 3)" << endl;
        exit(EXIT_FAILURE);
    }

    vector<Node> nodes;
    vector<Edge> edges;
    map<pair<double, double>, int> nodeId;
    map<int, vector<int>> adjList;

    // readFromRoadDhaka(file, nodes); // sir, amar ekta old function silo

    readDhakaStreet(argv[1], nodes, edges, nodeId, adjList);

    readDhakaRoute(argv[2], "metro", nodes, edges, nodeId, adjList);
    readDhakaRoute(argv[3], "bus", nodes, edges, nodeId, adjList);
    readDhakaRoute(argv[4], "bus", nodes, edges, nodeId, adjList);

    // Print the graph information (for verification)
    cout << "Number of nodes: " << nodes.size() << endl;
    cout << "Number of edges: " << edges.size() << endl;

    // for (const auto &edge : edges)
    // {
    //     cout << "Type: " << edge.type << ", Source: (" << nodes[edge.source].lat << ", "
    //          << nodes[edge.source].lng << ")";

    //     if (edge.type == "metro")
    //     {
    //         cout << " [" << edge.sourceName << "]"; // Only print station names for Metro
    //         cout << ", Destination: (" << nodes[edge.destination].lat << ", " << nodes[edge.destination].lng
    //              << ") [" << edge.destinationName << "]";
    //     }
    //     else
    //     {
    //         cout << ", Destination: (" << nodes[edge.destination].lat << ", " << nodes[edge.destination].lng
    //              << ")";
    //     }

    //     cout << ", Distance: " << edge.distance << " km, Altitude: " << edge.altitude << endl;
    // }

    double sourceLat, sourceLon, destLat, destLon;

    // sourceLat = 90.439973;
    // sourceLon = 23.741813;

    // destLat = 90.439341;
    // destLon = 23.742218;

    cout << "Enter the Source Latitude, Longitude, Destination Latitude, Longitude" << endl;
    cin >> sourceLat >> sourceLon >> destLat >> destLon;

    int srcNode = addNode(sourceLat, sourceLon, nodeId, nodes);
    int destNode = addNode(destLat, destLon, nodeId, nodes);

    // use dijkstra's algorithm to find the shortest path
    auto [path, totalDistance] = dijkstra(srcNode, destNode, nodes, edges, adjList);

    printf("Source: (%lf, %lf)\n", nodes[srcNode].lat, nodes[srcNode].lng);
    printf("Destination: (%lf, %lf)\n", nodes[destNode].lat, nodes[destNode].lng);

    cout << "\nShortest Path from (" << nodes[srcNode].lat << ", " << nodes[srcNode].lng << ") to ("
         << nodes[destNode].lat << ", " << nodes[destNode].lng << "):" << endl;
    bool hasPath = true;
    if (path.empty())
    {
        cout << "No path exists." << endl;
        hasPath = false;
    }
    else
    {
        cout << "Path: " << endl;
        for (size_t i = 0; i < path.size(); ++i)
        {
            if (i < path.size() - 1)
            {
                // Find Edge type and cost
                string edgeType = "None";
                string sourceName = "";
                string destinationName = "";
                double cost = 0.0;

                for (const auto &edge : edges)
                {
                    if (edge.source == path[i] && edge.destination == path[i + 1])
                    {
                        edgeType = edge.type;
                        sourceName = edge.sourceName;
                        destinationName = edge.destinationName;
                        cost = edge.distance * 2; // Assuming cost is distance * 2 for simplicity
                        break;
                    }
                }

                if (edgeType == "metro")
                {
                    cout << "Cost: ৳" << fixed << setprecision(2) << cost << ": Ride Metro from "
                         << fixed << setprecision(6)
                         << sourceName << " (" << nodes[path[i]].lat << ", " << nodes[path[i]].lng << ") to "
                         << destinationName << " (" << nodes[path[i + 1]].lat << ", " << nodes[path[i + 1]].lng << ")" << endl;
                }
                else if (edgeType == "road" || edgeType == "bus")
                {
                    cout << "Cost: ৳" << fixed << setprecision(2) << cost << ": Ride Bus from ("
                         << fixed << setprecision(6)
                         << nodes[path[i]].lat << ", " << nodes[path[i]].lng << ") to ("
                         << nodes[path[i + 1]].lat << ", " << nodes[path[i + 1]].lng << ")" << endl;
                }
                else
                {
                    cout << "Cost: ৳0.00: Walk from ("
                         << nodes[path[i]].lat << ", " << nodes[path[i]].lng << ") to ("
                         << nodes[path[i + 1]].lat << ", " << nodes[path[i + 1]].lng << ")" << endl;
                }
            }
        }
        cout << endl;

        cout << "Total Distance: " << totalDistance << " km" << endl;
    }

    if (hasPath){
        cout << "Generating KML file..." << endl;
        createKMLFile(nodes, edges, path);
        cout << "KML file generated successfully." << endl;
    }

    return 0;
}

// input 90.439973 23.741813 90.439341 23.742218
// input 90.344319 23.813352 90.342514 23.810823