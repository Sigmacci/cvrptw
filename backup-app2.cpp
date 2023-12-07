#include <time.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#define EXECUTION_TIME 60.0
#define RUN_FOREVER false               // set true if you want to run forever, but remember to set ITERATIONS_OF_GRASP to 0
#define ITERATIONS_OF_GRASP 280         // best tested value ??
#define RAND_PERCENTAGE 60              //  ,0-100 rand of greedy
#define ITERATIONS_OF_POINT_SEARCH 500  //

#define NUMBER_OF_ORIGINAL_PARENTS 10  //
#define NUMBER_OF_OFFSPRINGS 10        //
#define ITERATIONS_OF_OFFSPRINGS 100   //
#define GROWTH_FACTOR 1.1              //
#define PERCENTAGE_OF_MUTATION 40      // 0-100
#define PERCENTAGE_OF_CROSSOVER 10     // 0-100
#define MAX_ITERATIONS_WITHOUT_IMPROVEMENT 20

using namespace std;

typedef struct vehicle {
    int capacity;  // capacity of vehicle left
    double time;   // spent time
    int x_cord;
    int y_cord;

    vehicle(int capacity, int x_cord, int y_cord) {
        this->capacity = capacity;
        this->time = 0;
        this->x_cord = x_cord;
        this->y_cord = y_cord;
    };
} vehicle;
class Transport {
   public:
    int vehicle_cap;           // capacity of vehicle
    int x_cord_of_dispatcher;  // x coordinate of dispatcher
    int y_cord_of_dispatcher;  // y coordinate of dispatcher
    vector<vehicle> vehicles;  // storage of all vehicles

    Transport() {
        this->vehicle_cap = 0;
        this->x_cord_of_dispatcher = 0;
        this->y_cord_of_dispatcher = 0;
    }
    /// @brief create new vehicle and add it to vehicles vector
    /// @return index of new vehicle
    int dispatchNewVehicle() {
        vehicle new_vehicle(this->vehicle_cap, this->x_cord_of_dispatcher, this->y_cord_of_dispatcher);
        this->vehicles.push_back(new_vehicle);
        return this->vehicles.size() - 1;
    }
};

typedef struct customer {
    int id;
    int x_cord;
    int y_cord;
    int demand;
    int time_window_start;
    int time_window_end;
    int service_time;

    customer(int id, int x_cord, int y_cord, int demand, int time_window_start, int time_window_end, int service_time) {
        this->id = id;
        this->x_cord = x_cord;
        this->y_cord = y_cord;
        this->demand = demand;
        this->time_window_start = time_window_start;
        this->time_window_end = time_window_end;
        this->service_time = service_time;
    };
} customer;
class Customers {
   public:
    vector<customer> customers;  // storage of all customers and point 0
    Customers(){};
    /// @brief create new customer and add it to customers vector
    /// @return index of new customer
    int addCustomer(int id, int x_cord, int y_cord, int demand, int time_window_start, int time_window_end, int service_time) {
        customer new_customer(id, x_cord, y_cord, demand, time_window_start, time_window_end, service_time);
        this->customers.push_back(new_customer);
        return this->customers.size() - 1;
    }
};

class CodeExecutionCutoffTimer {
   public:
    chrono::time_point<chrono::system_clock> start;
    double targetTimeInSeconds;
    CodeExecutionCutoffTimer(double targetTimeInSeconds) {
        this->targetTimeInSeconds = targetTimeInSeconds;
        this->start = chrono::system_clock::now();
    }
    bool isTimeUp() {
        chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds = end - start;
        // return false;
        return elapsed_seconds.count() > targetTimeInSeconds;
    }
    double timeElapsed() {
        chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds = end - start;
        return elapsed_seconds.count();
    }
};

bool read_data_from_file(string path, Transport &transport, Customers &customers);
double time_to_arrive(int x1, int x2, int y1, int y2);
void saveToFile(vector<vector<int>> &solution, double cost);
void chechIfDone(CodeExecutionCutoffTimer &timer);
bool valid(Customers &customers, Transport &transport);
double cost(vector<vector<int>> &solution, Customers &customers, double **time_matrix);
void startFunction(Customers &customers, Transport &transport);
vector<vector<int>> greedy_randomized(Transport &transport, Customers &customers, double **time_matrix, double *best);

CodeExecutionCutoffTimer timer(EXECUTION_TIME);

// int ttt = 3;  //+2
// void tmp(int n, double time) {
//     ifstream file("solution.txt", ios::in);
//     ofstream file2("zzzout.txt", ios::app);
//     file2 << n << " " << file.rdbuf() << " " << time << endl;
//     file.close();
//     file2.close();
// }

int main(int argc, char *argv[]) {
    // check if file name is given
    if (argc != 2) {
        cout << argv[0] << " <file>" << endl;
        return 1;
    }

    // define variables for storing data
    Transport transport;
    Customers customers;

    // read data from file
    if (!read_data_from_file(argv[1], transport, customers)) {
        cout << "Error in reading data from file" << endl;
        return 1;
    }

    // check if data is valid
    if (!valid(customers, transport)) {
        ofstream file("solution.txt");
        file << -1 << " " << setprecision(5) << fixed << 0.0 << "\n";
        file.close();
        cout << "Solution does not exist" << endl;
        return 1;
    }

    startFunction(customers, transport);

    // for (int n = 27; n <= 1002; n += 25) {
    //     ttt = n;
    //     // define variables for storing data
    //     Transport transport;
    //     Customers customers;
    //     // read data from file
    //     if (!read_data_from_file(argv[1], transport, customers)) {
    //         cout << "Error in reading data from file" << endl;
    //         return 1;
    //     }

    //     if (!valid(customers, transport)) {
    //         ofstream file("solution.txt");
    //         file << -1 << " " << setprecision(5) << fixed << 0.0 << "\n";
    //         file.close();
    //         cout << "Solution does not exist" << endl;
    //         return 1;
    //     }
    //     cout << n << endl;
    //     chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();
    //     grasp(transport, customers);
    //     chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
    //     chrono::duration<double> elapsed_seconds = end - start;
    //     tmp(n - 2, elapsed_seconds.count());
    // }

    cout << "Solution in file solution.txt" << endl;
    return 0;
}

void chechIfDone(CodeExecutionCutoffTimer &t) {
    if (t.isTimeUp()) {
        cout << "Time is up" << endl
             << "Solution in file solution.txt" << endl;
        exit(0);
    }
}

/// @brief Read neccecery data from file and store them in vectors
/// @param path name and path of file
/// @return true->success false->fail
bool read_data_from_file(string path, Transport &transport, Customers &customers) {
    // open file
    ifstream file(path);
    if (!file.is_open()) {
        return false;
    }
    string line;
    getline(file, line);  // name
    getline(file, line);  // empty
    getline(file, line);  // VEHICLE
    getline(file, line);  // HEADERS
    getline(file, line);  // vehicle data
    istringstream iss(line);
    int tmp;
    iss >> tmp >> transport.vehicle_cap;
    getline(file, line);           // empty
    getline(file, line);           // CUSTOMER
    getline(file, line);           // HEADERS
    while (getline(file, line)) {  // customer data
                                   // while (getline(file, line) && ttt--) {  // customer data
        if (line.length() < 7)
            continue;
        istringstream iss(line);
        int id, x_cord, y_cord, demand, time_window_start, time_window_end, service_time;
        iss >> id >> x_cord >> y_cord >> demand >> time_window_start >> time_window_end >> service_time;
        if (id == 0) {
            transport.x_cord_of_dispatcher = x_cord;
            transport.y_cord_of_dispatcher = y_cord;
        }
        customers.addCustomer(id, x_cord, y_cord, demand, time_window_start, time_window_end, service_time);
    }
    file.close();
    return true;
}
void saveToFile(vector<vector<int>> &solution, double cost) {
    ofstream file("solution.txt");
    file << solution.size() << " " << setprecision(5) << fixed << cost << "\n";
    // file << solution.size() << " " << setprecision(5) << fixed << cost;
    for (int i = 0; i < solution.size(); i++) {
        for (int j = 0; j < solution[i].size(); j++) {
            file << solution[i][j] << " ";
        }
        file << "\n";
    }
    file.close();
}

bool valid(Customers &customers, Transport &transport) {
    for (customer customer : customers.customers) {
        if (customer.time_window_start > customer.time_window_end)
            return false;
        if (customer.time_window_start < 0)
            return false;
        // 0 -> dispatcher
        if (customer.time_window_start + customer.service_time >= customers.customers[0].time_window_end)
            return false;
        if (customer.service_time < 0)
            return false;
        if (customer.demand < 0)
            return false;
        if (customer.demand > transport.vehicle_cap)
            return false;
    }
    return true;
}

/// @brief Calculate time from vehicle to traget customer
/// @return time neccessary to arrive to customer
double time_to_arrive(int x1, int x2, int y1, int y2) {
    return sqrt(pow((double)x1 - (double)x2, 2) + pow((double)y1 - (double)y2, 2));
}

bool isRouteValid(vector<int> &route, int veh_cap, double **time_matrix, Customers &customers) {
    double time = 0.0;
    int cap = 0;
    for (int n = 0; n < route.size(); n++) {
        time += time_matrix[route[n]][route[n - 1]];
        if (time > (double)customers.customers[route[n]].time_window_end) {
            return false;
        }
        if (time <= (double)customers.customers[route[n]].time_window_start) {
            time = (double)customers.customers[route[n]].time_window_start;
        }
        time += (double)customers.customers[route[n]].service_time;
        cap += customers.customers[route[n]].demand;
        if (cap > veh_cap) {
            return false;
        }
    }
    if (time + time_matrix[route[route.size() - 1]][0] > (double)customers.customers[0].time_window_end) {
        return false;
    }
    return true;
}
double distance(int x1, int y1, int x2, int y2) {
    return sqrt(pow((double)x1 - (double)x2, 2) + pow((double)y1 - (double)y2, 2));
}
double cost(vector<vector<int>> &solution, Customers &customers, double **time_matrix) {
    double total_time = 0.0;
    for (int veh = 0; veh < solution.size(); veh++) {
        double veh_time = 0.0;
        int prev_point = 0;
        for (int point = 0; point < solution[veh].size(); point++) {
            double t = time_matrix[prev_point][solution[veh][point]];
            // double t = distance(customers.customers[solution[veh][point]].x_cord, customers.customers[solution[veh][point]].y_cord, customers.customers[prev_point].x_cord, customers.customers[prev_point].y_cord);
            if (veh_time + t > (double)customers.customers[solution[veh][point]].time_window_end) {
                return -1;
            }
            if (veh_time + t <= (double)customers.customers[solution[veh][point]].time_window_start) {
                veh_time = (double)customers.customers[solution[veh][point]].time_window_start;
            } else {
                veh_time += t;
            }
            veh_time += (double)customers.customers[solution[veh][point]].service_time;
            prev_point = solution[veh][point];
        }
        veh_time += time_matrix[prev_point][0];
        // veh_time += distance(customers.customers[0].x_cord, customers.customers[0].y_cord, customers.customers[prev_point].x_cord, customers.customers[prev_point].y_cord);
        if (veh_time > (double)customers.customers[0].time_window_end) {
            return -1;
        }
        total_time += veh_time;
    }
    // cout << "total: " << total_time << endl;
    return total_time;
}

bool is_acceptable(vector<vector<int>> &individual, int number_of_vehicles, int number_of_customers) {
    int check[number_of_customers] = {0};
    for (int i = 0; i < number_of_vehicles; i++) {
        for (int j = 0; j < number_of_customers; j++) {
            if (check[i] + individual[i][j] > 1) {
                individual[i][j] = 0;
            }
            check[i] += individual[i][j];
        }
        if (check[i] != 1) {
            return false;
        }
    }
    return true;
}

typedef struct {
    vector<vector<int>> routes;
    double cost;
} pop_result;

// TODO
bool isEmpty(pop_result &result) {
    for (int i = 0; i < result.routes.size(); i++)
        if (result.routes[i].empty())
            return true;
    return false;
}
void listRoutes(pop_result &result) {
    cout << "RESULT: " << endl;
    for (int i = 0; i < result.routes.size(); i++) {
        cout << i << ": " << result.routes[i].empty() << "  ";
        for (int j = 0; j < result.routes[i].size(); j++) {
            cout << result.routes[i][j] << " ";
        }
        cout << endl;
    }
}

bool check_if_violates_time_windows(vector<int> &route, Customers &customers, double **time_matrix) {
    double time = 0.0;
    int prev_point = 0;
    for (int i = 0; i < route.size(); i++) {
        time += time_matrix[prev_point][route[i]];
        if (time > (double)customers.customers[route[i]].time_window_end) {
            return true;
        }
        if (time <= (double)customers.customers[route[i]].time_window_start) {
            time = (double)customers.customers[route[i]].time_window_start;
        }
        time += (double)customers.customers[route[i]].service_time;
        prev_point = route[i];
    }
    time += time_matrix[prev_point][0];
    if (time > (double)customers.customers[0].time_window_end) {
        return true;
    }
    return false;
}
vector<int> try_new_routes(vector<int> &route, Customers &customers, double **time_matrix) {
    for (int i = 0; i < route.size(); i++) {
        for (int j = 0; j < route.size(); j++) {
            if (i == j) {
                continue;
            }
            vector<int> new_route;
            new_route.assign(route.begin(), route.end());
            int tmp = new_route[i];
            new_route.erase(new_route.begin() + i);
            new_route.insert(new_route.begin() + j, tmp);
            if (!check_if_violates_time_windows(new_route, customers, time_matrix)) {
                return new_route;
            }
        }
    }
    chechIfDone(timer);
    return vector<int>();
}
vector<vector<int>> fix_routes(vector<int> &route, Customers &customers, double **time_matrix, Transport &transport) {
    Transport tmp_transport;
    tmp_transport.vehicle_cap = transport.vehicle_cap;
    tmp_transport.x_cord_of_dispatcher = transport.x_cord_of_dispatcher;
    tmp_transport.y_cord_of_dispatcher = transport.y_cord_of_dispatcher;
    Customers tmp_customers;
    tmp_customers.addCustomer(0, customers.customers[0].x_cord, customers.customers[0].y_cord, 0, customers.customers[0].time_window_start, customers.customers[0].time_window_end, 0);
    for (int n = 0; n < route.size(); n++) {
        tmp_customers.addCustomer(customers.customers[route[n]].id, customers.customers[route[n]].x_cord, customers.customers[route[n]].y_cord, customers.customers[route[n]].demand, customers.customers[route[n]].time_window_start, customers.customers[route[n]].time_window_end, customers.customers[route[n]].service_time);
    }
    double best;
    vector<vector<int>> response = greedy_randomized(tmp_transport, tmp_customers, time_matrix, &best);
    pop_result tmp;
    tmp.routes.clear();
    tmp.routes.assign(response.begin(), response.end());
    return response;
}

vector<vector<int>> greedy_randomized(Transport &transport, Customers &customers, double **time_matrix, double *best) {
    vector<vector<int>> result;
    vector<int> route;
    vector<double> coefficients;
    vector<customer> customers_copy;
    customers_copy.assign(customers.customers.begin() + 1, customers.customers.end());
    int j = transport.dispatchNewVehicle(), r = 0;
    while (!customers_copy.empty()) {
        int prev_customer_id = 0;
        if (!route.empty()) {
            prev_customer_id = route.back();
        }
        double max = 0;
        int maxk = -1;
        int k = 0;
        for (; k < customers_copy.size(); k++) {
            coefficients.push_back((double)customers_copy[k].demand / time_matrix[prev_customer_id][customers.customers[k].id]);
            if (max < coefficients.back()) {
                max = coefficients.back();
                maxk = k;
            };
        }
        int iter = 0;
        while (++iter < ITERATIONS_OF_POINT_SEARCH) {
            if (rand() % 101 <= RAND_PERCENTAGE && customers_copy.size() > 1)
                r = (rand() % (customers_copy.size() - 1)) + 1;
            else
                r = maxk;

            double tmpT = time_matrix[prev_customer_id][customers_copy[r].id] + transport.vehicles[j].time;
            if (tmpT <= (double)customers_copy[r].time_window_end && transport.vehicles[j].capacity >= customers_copy[r].demand && tmpT + (double)customers_copy[r].service_time + time_matrix[customers_copy[r].id][0] <= (double)customers.customers[0].time_window_end) {
                transport.vehicles[j].time = tmpT < (double)customers_copy[r].time_window_start ? (double)customers_copy[r].time_window_start : tmpT;
                transport.vehicles[j].time += (double)customers_copy[r].service_time;
                transport.vehicles[j].x_cord = customers_copy[r].x_cord;
                transport.vehicles[j].y_cord = customers_copy[r].y_cord;
                transport.vehicles[j].capacity -= customers_copy[r].demand;
                route.push_back(customers_copy[r].id);
                customers_copy.erase(customers_copy.begin() + r);
                coefficients.clear();
                iter = ITERATIONS_OF_POINT_SEARCH + 1;
            }
            chechIfDone(timer);
        }
        if (iter < ITERATIONS_OF_POINT_SEARCH + 1) {
            result.push_back(route);
            route.clear();
            j = transport.dispatchNewVehicle();
        }
    }
    result.push_back(route);
    *best = cost(result, customers, time_matrix);
    return result;
}

void fix_crossed_child(pop_result &child, Customers &customers, double **time_matrix, Transport &transport) {
    bool points[customers.customers.size() - 1] = {false};
    int tmp = 0;
    for (int n = child.routes.size() - 1; n >= 0; n--) {
        for (int m = child.routes[n].size() - 1; m >= 0; m--) {
            int tmp = child.routes[n][m];
            if (points[child.routes[n][m] - 1]) {
                child.routes[n].erase(child.routes[n].begin() + m);
                if (child.routes[n].empty()) {
                    child.routes.erase(child.routes.begin() + n);
                }
            }
            points[tmp - 1] = true;
        }
    }
    vector<int> missing_points;
    for (int n = 0; n < customers.customers.size() - 1; n++)
        if (!points[n])
            missing_points.push_back(n + 1);
    if (missing_points.empty()) return;
    vector<vector<int>> result = fix_routes(missing_points, customers, time_matrix, transport);

    for (int i = 0; i < result.size(); i++)
        child.routes.push_back(result[i]);
}
void mutate(pop_result &child, Customers &customers, double **time_matrix) {
    for (int k = 0; k < (PERCENTAGE_OF_MUTATION * child.routes.size()) / 100; k++) {
        // choose random car to mutate
        int car = rand() % child.routes.size();
        if (child.routes[car].size() >= 2) {
            vector<int> mutatedCar;
            mutatedCar.assign(child.routes[car].begin(), child.routes[car].end());

            int n = rand() % mutatedCar.size();
            int m = rand() % mutatedCar.size();
            while (n == m) {
                m = rand() % mutatedCar.size();
            }
            int tmp = mutatedCar[n];
            mutatedCar[n] = mutatedCar[m];
            mutatedCar[m] = tmp;
            // validate and fix the mutated car
            if (check_if_violates_time_windows(mutatedCar, customers, time_matrix)) {
                mutatedCar.clear();
                mutatedCar = try_new_routes(child.routes[car], customers, time_matrix);
            }
            if (mutatedCar.empty())
                continue;
            // save the mutated car
            child.routes[car].clear();
            child.routes[car].assign(mutatedCar.begin(), mutatedCar.end());
        }
    }
}
bool save_child(pop_result &child, pop_result &bestResult, vector<pop_result> &population, int targetSize, Customers &customers, double **time_matrix) {
    child.cost = cost(child.routes, customers, time_matrix);
    cout << child.cost << endl;
    if (child.cost != -1) {
        bool found = false;
        for (int i = 0; i < population.size(); i++) {
            if (population[i].cost > child.cost) {
                population.insert(population.begin() + i, child);
                found = true;
                break;
            }
        }
        if (!found)
            population.push_back(child);
        if (population.size() >= targetSize)
            population.pop_back();
        if (child.cost < bestResult.cost) {
            bestResult.routes.clear();
            bestResult.routes.assign(child.routes.begin(), child.routes.end());
            bestResult.cost = child.cost;
            saveToFile(bestResult.routes, bestResult.cost);
            return true;
        }
    }
    return false;
}

pop_result genetic_algorithm(Transport &transport, Customers &customers, double **time_matrix) {
    vector<pop_result> population;
    // start by generating one parent
    pop_result bestResult;
    bestResult.cost = numeric_limits<double>::max();
    for (int i = 0; i < NUMBER_OF_ORIGINAL_PARENTS; i++) {
        pop_result parent;
        parent.routes = greedy_randomized(transport, customers, time_matrix, &parent.cost);
        if (parent.cost < bestResult.cost) {
            bestResult.cost = parent.cost;
            bestResult.routes.clear();
            bestResult.routes.assign(parent.routes.begin(), parent.routes.end());
            saveToFile(bestResult.routes, bestResult.cost);
        }
        population.push_back(parent);
    }

    int iter = 0;
    int iterations_without_improvement = 0;
    while (iter++ < ITERATIONS_OF_OFFSPRINGS && ++iterations_without_improvement <= MAX_ITERATIONS_WITHOUT_IMPROVEMENT) {
        int targetSize = (GROWTH_FACTOR * population.size());
        cout << "iteration: " << iter << " pop size: " << population.size() << " target size: " << targetSize << endl;
        vector<pop_result> new_population;
        // generate x offsprings from current population
        for (int parent = 0; parent < population.size(); parent++) {
            save_child(population[parent], bestResult, new_population, targetSize, customers, time_matrix);
            for (int i = 0; i < NUMBER_OF_OFFSPRINGS; i++) {
                // select random father
                int father = rand() % population.size();
                while (father == parent)
                    father = rand() % population.size();

                pop_result mother_child, father_child;
                mother_child.routes.assign(population[parent].routes.begin(), population[parent].routes.end());
                father_child.routes.assign(population[father].routes.begin(), population[father].routes.end());

                // crossover
                int minGeneticSize = min(mother_child.routes.size(), father_child.routes.size());
                for (int p = 0; p < (PERCENTAGE_OF_CROSSOVER * minGeneticSize) / 100; p++) {
                    int car1 = rand() % mother_child.routes.size();
                    int car2 = rand() % father_child.routes.size();
                    vector<int> newCar1, newCar2;
                    int smallestSize = min(mother_child.routes[car1].size(), father_child.routes[car2].size());

                    int crossover_point_left = rand() % smallestSize;                                                         //<0;size)
                    int crossover_point_right = crossover_point_left + (rand() % (smallestSize - crossover_point_left)) + 1;  //(left;size>
                    // crossover
                    for (int n = 0; n < crossover_point_left; n++) {
                        newCar1.push_back(mother_child.routes[car1][n]);
                        newCar2.push_back(father_child.routes[car2][n]);
                    }
                    for (int n = crossover_point_left; n < crossover_point_right; n++) {
                        newCar1.push_back(father_child.routes[car2][n]);
                        newCar2.push_back(mother_child.routes[car1][n]);
                    }
                    for (int n = crossover_point_right; n < mother_child.routes[car1].size(); n++)
                        newCar1.push_back(mother_child.routes[car1][n]);

                    for (int n = crossover_point_right; n < father_child.routes[car2].size(); n++)
                        newCar2.push_back(father_child.routes[car2][n]);

                    // fix the child if it is not valid
                    mother_child.routes[car1].clear();
                    if (!isRouteValid(newCar1, transport.vehicle_cap, time_matrix, customers)) {
                        vector<int> tmp = try_new_routes(newCar1, customers, time_matrix);
                        vector<vector<int>> fixed;
                        if (tmp.size() == 0) {
                            fixed = fix_routes(newCar1, customers, time_matrix, transport);
                            mother_child.routes[car1].assign(fixed[0].begin(), fixed[0].end());
                            for (int i = 1; i < fixed.size(); i++)
                                mother_child.routes.push_back(fixed[i]);
                        } else {
                            mother_child.routes[car1].assign(tmp.begin(), tmp.end());
                        }
                    } else {
                        mother_child.routes[car1].assign(newCar1.begin(), newCar1.end());
                    }

                    father_child.routes[car2].clear();
                    if (!isRouteValid(newCar2, transport.vehicle_cap, time_matrix, customers)) {
                        vector<int> tmp = try_new_routes(newCar2, customers, time_matrix);
                        vector<vector<int>> fixed;
                        if (tmp.size() == 0) {
                            fixed = fix_routes(newCar1, customers, time_matrix, transport);
                            father_child.routes[car2].assign(fixed[0].begin(), fixed[0].end());
                            for (int i = 1; i < fixed.size(); i++)
                                father_child.routes.push_back(fixed[i]);
                        } else {
                            father_child.routes[car2].assign(tmp.begin(), tmp.end());
                        }
                    } else {
                        father_child.routes[car2].assign(newCar2.begin(), newCar2.end());
                    }
                    // final fix
                    fix_crossed_child(mother_child, customers, time_matrix, transport);
                    fix_crossed_child(father_child, customers, time_matrix, transport);
                }
                // mutation

                mutate(mother_child, customers, time_matrix);
                mutate(father_child, customers, time_matrix);

                // select the best child
                save_child(mother_child, bestResult, new_population, targetSize, customers, time_matrix);
                save_child(father_child, bestResult, new_population, targetSize, customers, time_matrix);

                chechIfDone(timer);
            }
        }
        population.clear();
        population.assign(new_population.begin(), new_population.end());
    }
    return bestResult;
}

void startFunction(Customers &customers, Transport &transport) {
    srand(time(NULL));
    // generate time matrix
    double **time_matrix = new double *[customers.customers.size()];
    for (int i = 0; i < customers.customers.size(); i++) {
        time_matrix[i] = new double[customers.customers.size()];
        for (int j = 0; j < i + 1; j++) {
            time_matrix[i][j] = time_to_arrive(customers.customers[i].x_cord, customers.customers[j].x_cord, customers.customers[i].y_cord, customers.customers[j].y_cord);
            time_matrix[j][i] = time_matrix[i][j];
        }
        chechIfDone(timer);
    }
    // start copulating
    genetic_algorithm(transport, customers, time_matrix);
}