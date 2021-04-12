#include<iostream>
#include <fstream>
#include <string>
#include<set>
#include<list>
#include<map>
#include<cstdlib>
#include<ctime>

typedef struct vertex {
    int index;
    bool is_visited = false;
    std::set<struct vertex *> neighbors;
} vertex;

//Function Decleration:
std::set<vertex *> create_graph(int, double);

void connect_Vertices(vertex *, vertex *);

bool chance_of_connecting(double p);

int is_isolated(std::set<vertex *> &graph);

std::map<int, int> distance_from_vertex(vertex *s);

int diameter(std::set<vertex *> &graph);

int connectivity(std::set<vertex *> &graph);

std::set<vertex *> create_graph(int n, double p) {
    int count = 0;
    srand(time(NULL));
    std::set<vertex *> graph;
    for (int k = 0; k < n; k++) {
        graph.insert(new vertex{k});
        // std::cout << "index number " << k << " has been created" << std::endl;
    }
    for (std::set<vertex *>::iterator i = graph.begin(); i != graph.end(); ++i) {
        for (std::set<vertex *>::iterator j = next(i); j != graph.end(); ++j) {
            if (chance_of_connecting(p)) {
                connect_Vertices(*i, *j);
                count++;
                //std::cout << (*i)->index << "--------------->" << (*j)->index << std::endl;
            }

        }
    }
    //std::cout << "the number of edges is " << count << std::endl;
    return graph;
}

void connect_Vertices(vertex *a, vertex *b) {
    a->neighbors.insert(b);
    b->neighbors.insert(a);
}

bool chance_of_connecting(double p) {
    return (RAND_MAX * p >= rand());
}

int is_isolated(std::set<vertex *> &graph) {
    for (std::set<vertex *>::iterator i = graph.begin(); i != graph.end(); ++i) {
        if ((*i)->neighbors.size() == 0)
            return 1;
    }
    return 0;
}

std::map<int, int> distance_from_vertex(vertex *s) {
    bool flag = !(s->is_visited);
    s->is_visited = flag;
    std::list<vertex *> queue;
    std::map<int, int> distance_mapping; //vertex index, distance
    distance_mapping[s->index] = 0;
    queue.push_back(s);
    while (!(queue.empty())) {
        s = queue.front();
        queue.pop_front();
        for (std::set<vertex *>::iterator i = s->neighbors.begin(); i != s->neighbors.end(); ++i) {
            if ((*i)->is_visited != flag) {
                queue.push_back(*i);
                (*i)->is_visited = flag;
                distance_mapping[(*i)->index] = distance_mapping[s->index] + 1;
            }
        }
    }
    return distance_mapping;
}

int diameter(std::set<vertex *> &graph) {
    int max = 0;
    for (std::set<vertex *>::iterator i = graph.begin(); i != graph.end(); i++) {
        std::map<int, int> distance_map = distance_from_vertex(*i);
        for (std::map<int, int>::iterator j = distance_map.begin(); j != distance_map.end(); ++j) {
            if (j->second > max)
                max = j->second;
        }
    }
    return max;
}

int connectivity(std::set<vertex *> &graph) {
    if ((distance_from_vertex(*graph.begin()).size()) == graph.size())
        return 1;
    return 0;
}

void predictions_to_csv() {
    int v = 1000, presicion = 500, connectivity_counter = 0, diameter_counter = 0, is_isolated_counter = 0, progress;
    double treshhold1 = (log(v)) / v, treshhold2 = sqrt((2 * log(v)) / v);
    std::ofstream file("graph probability.csv");
    for (int condition = 1; condition <= 3; condition++) {
        progress = 0;
        std::cout << "Condition " << condition << std::endl;
        if (condition == 1) {
            file << "Connectivity,";
            file << "p,";
            for (double p = ((treshhold1 / 5) - (treshhold1 / 10)); p < 2 * treshhold1; p += treshhold1 / 5) {
                file << std::to_string(p) << ",";
                connectivity_counter = 0;
                for (int i = 0; i < presicion; i++) {
                    std::set<vertex *> graph = create_graph(v, p);
                    if ((p < treshhold1) ^ connectivity(graph)) {
                        connectivity_counter++;
                    }
                    for (std::set<vertex *>::iterator j = graph.begin(); j != graph.end(); j++) {
                        delete *j;
                    }
                    if (!(i % (50)))
                        std::cout << "\r" << ++progress << "%";
                }
                file << (double) connectivity_counter / presicion << ",";
            }
        }
        if (condition == 2) {

            file << "Diameter,";
            file << "p,";
            for (double p = ((treshhold2 / 5) - (treshhold2 / 10)); p < 2 * treshhold2; p += treshhold2 / 5) {
                file << std::to_string(p) << ",";
                diameter_counter = 0;
                for (int i = 0; i < presicion; i++) {
                    std::set<vertex *> graph = create_graph(v, p);
                    if(connectivity(graph)){
                        if ((p < treshhold2) ^ (diameter(graph) == 2)) {
                            diameter_counter++;
                        }
                    }
                    for (std::set<vertex *>::iterator j = graph.begin(); j != graph.end(); j++) {
                        delete *j;
                    }
                    if (!(i % (50)))
                        std::cout << "\r" << ++progress << "%";
                }
                file << (double) diameter_counter / presicion << ",";
            }
        }
        if (condition == 3) {

            file << "Is Isolated,";
            file << "p,";
            for (double p = ((treshhold1 / 5) - (treshhold1 / 10)); p < 2 * treshhold1; p += treshhold1 / 5) {
                file << std::to_string(p) << ",";
                is_isolated_counter = 0;
                for (int i = 0; i < presicion; i++) {
                    std::set<vertex *> graph = create_graph(v, p);
                    if ((p > treshhold1) ^ (is_isolated(graph))) {
                        is_isolated_counter++;
                    }
                    for (std::set<vertex *>::iterator j = graph.begin(); j != graph.end(); j++) {
                        delete *j;
                    }
                    if (!(i % (50)))
                        std::cout << "\r" << ++progress << "%";
                }
                file << (double) is_isolated_counter / presicion << ",";
            }
        }
    }
    file.close();
}

int main() {
    predictions_to_csv();

/*
    std::set<vertex *> graph = create_graph(7, 0.5);
    if (is_isolated(graph) == 1)
        std::cout << "Graph has an isolated vertex" << std::endl;
    else
        std::cout << "Graph has NO isolated vertex" << std::endl;

    std::cout << "the diameter of the graph is: " << diameter(graph) << std::endl;

    if (connectivity(graph))
        std::cout << "the graph is connected" << std::endl;
    else {
        std::cout << "the graph is NOT connected" << std::endl;
    }
    */
}