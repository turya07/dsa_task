#include <bits/stdc++.h>
#define INT unsigned int

#define INF (INT)1e8
using namespace std;

struct Edge
{
    std::string edgeName;
    std::string nodeA;
    std::string nodeB;
    double weight;
};

int main(const int argc, const char *argv[])
{

    if (argc != 2)
    {
        cerr << "Usage: " << "./main" << " <data-set>" << endl;
        exit(EXIT_FAILURE);
    }
    ifstream file;
    file.open(argv[1], ios::out);
    if (!file.is_open())
    {
        cerr << "Failed to load data set" << endl;
        exit(EXIT_FAILURE);
        return 1;
    }
    int n = 0;
    string s, line;
    int tracking = 0;

    int i = 0;
    vector<Edge> edges;
    while (getline(file, line))
    {
        i++;
        
        char comma = ',';
        stringstream ss(line);
        string roadName;
        getline(ss, roadName, comma);

        float lon1, lat1, lon2, lat2;
        stack<string> weight;
        ss >> lon1 >> comma >> lat1 >> comma >> lon2 >> comma >> lat2 >> comma;

        while (ss.good())
        {
            string substr;
            getline(ss, substr, ',');
            weight.push(substr);
        }

        edges.push_back({roadName, to_string(lon1) + ", " + to_string(lat1), to_string(lon2) + ", " + to_string(lat2), double_t(stod(weight.top()))});

        // poping the last 2 element of the stack
        weight.pop();
        weight.pop();
        
        // optimizing edge and reaching to the last pari of LNG and LAT
        while (!weight.empty())
        {
            if (weight.size() == 2)
            {
                lat2 = stod(weight.top());
                weight.pop();
                lon2 = stod(weight.top());
                weight.pop();
            }
            break;
        }

        Edge *top;
        top = &edges.back();
        top->nodeB = to_string(lon2) + ", " + to_string(lat2); // updating the last pair of LNG and LAT
        cout << "Name: " << roadName << endl;
        cout << "Source: " << "(" << top->nodeA << ")" << endl;
        cout << "Destination: " << "(" << top->nodeB << ")" << endl;

        cout << "Distance in Meters: " << fixed << setprecision(2) << top->weight * 1000 << endl
             << endl;
    }
    cout << "Number of Data: " << i << endl
         << endl;
    return 0;

    vector<vector<INT>> graph;
    graph.resize(n, vector<INT>(n));

    for (INT i = 0; i < n; i++)
    {
        for (INT j = 0; j < n; j++)
        {
            cin >> graph[i][j];
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