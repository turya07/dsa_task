#include <bits/stdc++.h>
using namespace std;

#define INT unsigned int
#define INF (INT)1e8

typedef struct Node Node;
typedef struct Edge Edge;
struct Node
{
    double lat;
    double lng;
    string lat_lng;
    vector<Edge> edges;
};

struct Edge
{
    Node dest;
    double weight;
};
bool operator<(const Node &a, const Node &b)
{
    return a.lat < b.lat;
}

bool operator==(const Node &a, const Node &b)
{
    return a.lat_lng == b.lat_lng;
}
void readDataFromFile(ifstream &file, set<Node> &nodes)
{
    int counter = 0;
    string line;
    while (getline(file, line))
    {
        counter++;
        char comma = ',';
        stringstream ss(line);
        vector<Edge> edges;
        string roadName;
        getline(ss, roadName, comma);

        double lon1, lat1, lon2, lat2;
        stack<string> weight;
        ss >> lon1 >> comma >> lat1 >> comma;
        // eituku obdi ami proti line read korsi, then first node er jonno lat,lng 1 and 2nd node er jonno lat,lng 2 nisi.
        Node a, b, *temp;
        a.lat = lat1;
        a.lng = lon1;
        a.lat_lng = to_string(lat1) + "," + to_string(lon1);

        while (ss.good())
        {
            // ei loop er maddhome ekdom last obdi read krsi, last e ekta distance/weight dewa ase, oitake weight stack e push krsi ekdom finally.
            ss >> lon2 >> comma >> lat2 >> comma;
            weight.push(to_string(lon2));
            weight.push(to_string(lat2));
        }

        // cout << ": " << weight.top() << endl;
        double w = stod(weight.top());
        weight.pop();
        weight.pop();

        while (!weight.empty())
        {
            lat2 = stod(weight.top());
            weight.pop();
            lon2 = stod(weight.top());
            weight.pop();
            b.lat = lat2;
            b.lng = lon2;
            b.lat_lng = to_string(lat2) + "," + to_string(lon2);
            edges.push_back({b, w});
        }

        // check krtesi je first node ta er age amar SET of Node e ase kina
        if (nodes.find(a) != nodes.end())
        {
            a = *const_cast<Node *>(&(*nodes.find(a))); // jodi thake, temp = node.find(a) krsi.
        }
        for (auto e : edges)
        {
            a.edges.push_back(e); // temp er edges e `e node` push krsi.
        }

        nodes.insert(a); // finally temp er value ta insert krsi nodes set e. ekhane SET use koray automatic duplicate value thakbe na. But old destination node ta thakle new destination node ta soho new destination node add hobe edge e.
        edges.clear();
    }
}

// main function
int main(const int argc, const char *argv[])
{

    // arguement theke file nicchi, tai argc == 2 hote hobe
    if (argc != 2)
    {
        cerr << "Usage: " << "./main" << " <data-set>" << endl;
        exit(EXIT_FAILURE);
    }
    ifstream file;
    ofstream output;
    file.open(argv[1], ios::out);
    if (!file.is_open())
    {
        cerr << "Failed to load data set" << endl;
        exit(EXIT_FAILURE);
        return 1;
    }

    set<Node> nodes;

    readDataFromFile(file, nodes);
    // now write data in output terminal

    output.open("output.txt", ios::out);
    unsigned int tt = 0, te = 0;
    for (auto node : nodes)
    {
        output << node.lat_lng << ":";
        for (auto edge : node.edges)
        {
            output << "\t=>" << edge.dest.lat_lng;
        }
        output << "--> " + to_string(node.edges.front().weight) << endl;
    }
    output << "========++========++========++========++========++========" << endl;
    output.close();

    return 0;

    // Nicher code ta normal ekta dijkstra algorithm er graph

    int n = 0;
    vector<vector<int>> graph;
    graph.resize(n, vector<int>(n));

    for (INT i = 0; i < n; i++)
    {
        for (INT j = 0; j < n; j++)
        {
            double lat1, lat2, lon1, lon2;
            cin >> lat1 >> lon1 >> lat2 >> lon2;
        }
    }

    for (INT i = 0; i < n; i++)
    {
        cout << (i) << ": ";
        for (INT j = 0; j < n; j++)
        {
            if (graph[i][j] > 0)
            {
                cout << (j) << "(" << graph[i][j] << ") ";
            }
        }
        cout << endl;
    }

    vector<INT> dist(n, INF);
    vector<INT> parent(n, -1);
    vector<bool> visited(n, false);

    dist[0] = 0;
    for (INT i = 0; i < n; i++)
    {
        INT min_dist = INF;
        INT min_index = -1;
        for (INT j = 0; j < n; j++)
        {
            if (!visited[j] && dist[j] < min_dist)
            {
                min_dist = dist[j];
                min_index = j;
            }
        }

        visited[min_index] = true;

        for (INT j = 0; j < n; j++)
        {
            if (!visited[j] && graph[min_index][j] > 0 && dist[j] > dist[min_index] + graph[min_index][j])
            {
                dist[j] = dist[min_index] + graph[min_index][j];
                parent[j] = min_index;
            }
        }
    }

    // Print the shortest path
    for (INT i = 0; i < n; i++)
    {
        cout << "Shortest path from 1 to " << i + 1 << " is: ";
        INT j = i;
        while (j != -1)
        {
            cout << j + 1 << " ";
            j = parent[j];
        }
        cout << "with cost: " << dist[i] << endl;
    }

    return 0;
}