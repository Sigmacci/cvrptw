#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

typedef struct vehicle
{
    int capacity;
    int time;
    int x_cord;
    int y_cord;

    vehicle(int capacity)
    {
        this->capacity = capacity;
        this->time = 0;
        this->x_cord = 0;
        this->y_cord = 0;
    };
} vehicle;

typedef struct receiver
{
    int id;
    int x_cord;
    int y_cord;
    int demand;
    int time_window_start;
    int time_window_end;
    int service_time;

    receiver(int id, int x_cord, int y_cord, int demand, int time_window_start, int time_window_end, int service_time)
    {
        this->id = id;
        this->x_cord = x_cord;
        this->y_cord = y_cord;
        this->demand = demand;
        this->time_window_start = time_window_start;
        this->time_window_end = time_window_end;
        this->service_time = service_time;
    };
} receiver;

int time_to_arrive(vehicle v, receiver r);
void heapify(std::vector<receiver> receivers, std::vector<double> coefficient, int n, int i);
void heapsort(std::vector<receiver> receivers, std::vector<double> coefficient, int n);
std::vector<std::vector<int>> greedy_randomized(std::vector<vehicle> vehicles, std::vector<receiver> receivers);

int main(int argc, char *argv[])
{
    std::vector<vehicle> vehicles;
    std::vector<receiver> receivers;
    receiver home(0, 0, 0, 0, 0, 0, 0);
    receivers.push_back(home);
}

int time_to_arrive(vehicle v, receiver r)
{
    int t = static_cast<int>(sqrt(
        pow(v.x_cord - r.x_cord, 2) +
        pow(v.y_cord - r.y_cord, 2)));

    if (v.time + t < r.time_window_start)
        return r.time_window_start - v.time;

    return t;
}

void heapify(std::vector<receiver> receivers, std::vector<double> coefficient, int n, int i)
{
    int min = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    if (left < n && coefficient[left] < coefficient[min])
        min = left;
    if (right < n && coefficient[right] < coefficient[min])
        min = right;
    if (min != i)
    {
        std::swap(coefficient[i], coefficient[min]);
        std::swap(receivers[i], receivers[min]);
        heapify(receivers, coefficient, n, min);
    }
}

void heapsort(std::vector<receiver> receivers, std::vector<double> coefficient, int n)
{
    for (int i = n / 2 - 1; i >= 0; i--)
    {
        heapify(receivers, coefficient, n, i);
    }
    for (int i = 1; i < n; i++)
    {
        std::swap(coefficient[0], coefficient[n - i]);
        std::swap(receivers[0], receivers[n - i]);
        heapify(receivers, coefficient, n - i, 0);
    }
}

std::vector<std::vector<int>> greedy_randomized(std::vector<vehicle> vehicles, std::vector<receiver> receivers)
{
    std::vector<std::vector<int>> result;
    std::vector<int> route = {0};
    route.resize(receivers.size() - 1);
    std::vector<double> coefficients;
    std::vector<vehicle> vehicles_copy;
    std::vector<receiver> receivers_copy;
    std::copy(vehicles, vehicles_copy, vehicles.size());
    std::copy(receivers, receivers_copy, receivers.size());

    int i = 1, j = 0, r = 0;

    while (i < receivers_copy.size())
    {
        for (int k = 0; k < receivers_copy.size(); k++)
        {
            coefficients.push_back(receivers_copy[k].demand / time_to_arrive(vehicles_copy[j], receivers_copy[k]));
        }
        heapsort(receivers_copy, coefficients, receivers_copy.size());

        if (rand() % 2 == 0)
            r = (rand() % (receivers_copy.size() - 1)) + 1;
        else
            r = i;

        int t = time_to_arrive(vehicles_copy[j], receivers_copy[r]);
        if (vehicles_copy[j].capacity >= receivers_copy[r].demand &&
            vehicles_copy[j].time + t <= receivers_copy[r].time_window_end)
        {
            vehicles_copy[j].capacity -= receivers_copy[r].demand;
            vehicles_copy[j].time += t + receivers_copy[r].service_time;
            vehicles_copy[j].x_cord = receivers_copy[r].x_cord;
            vehicles_copy[j].y_cord = receivers_copy[r].y_cord;
            route[receivers_copy[r].id - 1] = 1;
            receivers_copy.erase(receivers_copy.begin() + i);
            coefficients.clear();
        }
        else
        {
            if (i == receivers_copy.size())
            {
                i = 0;
                result.push_back(route);
                if (++j == vehicles_copy.size())
                {
                    return {{}};
                }
                continue;
            }
        }
        if (r == i)
            i++;
    }
    return result;
}

int fitness(std::vector<std::vector<int>> individual, std::vector<receiver> receivers, int capacity)
{
    int total_fitness = 0;
    int temp_fitness = 0;
    vehicle dummy_vehicle(capacity);
    for (auto route : individual)
    {
        for (int i = 1; i < receivers.size(); i++)
        {
            if (route[i] == 1)
            {
                int t = time_to_arrive(dummy_vehicle, receivers[i]);
                if (t > receivers[i].time_window_end)
                {
                    return -1;
                }
                dummy_vehicle.capacity -= receivers[i].demand;
                dummy_vehicle.time += t + receivers[i].service_time;
                dummy_vehicle.x_cord = receivers[i].x_cord;
                dummy_vehicle.y_cord = receivers[i].y_cord;
                temp_fitness = dummy_vehicle.time;
            }
        }
        dummy_vehicle.capacity = capacity;
        dummy_vehicle.time = 0;
        dummy_vehicle.x_cord = 0;
        dummy_vehicle.y_cord = 0;
        total_fitness += temp_fitness;
    }
    return total_fitness;
}

bool is_acceptable(std::vector<std::vector<int>> individual, int number_of_vehicles, int number_of_receivers)
{
    int check[number_of_receivers] = {0};
    for (int i = 0; i < number_of_vehicles; i++)
    {
        for (int j = 0; j < number_of_receivers; j++)
        {
            if (check[i] + individual[i][j] > 1)
            {
               individual[i][j] = 0;
            }
            check[i] += individual[i][j];
        }
        if (check[i] != 1)
        {
            return false;
        }
    }
    return true;
}

int genetic_algorithm(std::vector<vehicle> vehicles, std::vector<receiver> receivers, int capacity, int number_of_individuals, float mutation_rate)
{
    int k = 0;
    std::vector<std::vector<std::vector<int>>> population;
    std::vector<int> fitnesses;
    int min = -1;
    for (int i = 0; i < number_of_individuals; i++)
    {
        population.push_back(greedy_randomized(vehicles, receivers));
        fitnesses.push_back(fitness(population[i], receivers, capacity));
        min = fitnesses[i] < min || min == -1 ? fitnesses[i] : min;
    }
    while (k < 100)
    {
        std::vector<std::vector<std::vector<int>>> new_population;
        for (int i = 0; i < number_of_individuals; i++)
        {
            fitnesses.clear();
            // parent selection
            int parent1 = rand() % number_of_individuals;
            int parent2 = rand() % number_of_individuals;
            while (parent1 == parent2)
            {
                parent2 = rand() % number_of_individuals;
            }

            int offspring1[vehicles.size()][receivers.size() - 1], offspring2[vehicles.size()][receivers.size() - 1];
            std::vector<std::vector<int>> final_offspring;
            final_offspring.resize(vehicles.size());

            // crossover
            for (int j = 0; j < vehicles.size(); j++)
            {
                final_offspring[j].resize(receivers.size() - 1);
                for (int l = 0; l < receivers.size() - 1; l++)
                {
                    if (population[parent1][j][l] == population[parent2][j][l])
                    {
                        offspring1[j][l] = population[parent1][j][l];
                        offspring2[j][l] = population[parent1][j][l];
                    }
                    else
                    {
                        if (rand() % 2 == 0)
                        {
                            offspring1[j][l] = population[parent1][j][l];
                            offspring2[j][l] = population[parent2][j][l];
                        }
                        else
                        {
                            offspring1[j][l] = population[parent2][j][l];
                            offspring2[j][l] = population[parent1][j][l];
                        }
                    }
                }
                int decimal1 = 0, decimal2 = 0;
                for (int l = receivers.size() - 2; l >= 0; l--)
                {
                    decimal1 += pow(2, receivers.size() - l - 2) * offspring1[j][l];
                    decimal2 += pow(2, receivers.size() - l - 2) * offspring2[j][l];
                }
                decimal1 = receivers.size() - decimal1 - 1;
                decimal2 = receivers.size() - decimal2 - 1;
                for (int l = receivers.size() - 2; l >= 0; l--)
                {
                    offspring1[j][l] = decimal1 % 2;
                    decimal1 /= 2;
                    offspring2[j][l] = decimal2 % 2;
                    decimal2 /= 2;
                    final_offspring[j][l] = (offspring1[j][l] & population[parent1][j][l]) | (offspring2[j][l] & population[parent2][j][l]);

                    //mutation
                    if (rand() % 100 < mutation_rate)
                    {
                        final_offspring[j][l] = ~final_offspring[j][l];
                    }
                }
                int t = fitness(final_offspring, receivers, capacity);
                if (is_acceptable(final_offspring, vehicles.size(), receivers.size() - 1) && t != -1)
                {
                    new_population[i].push_back(final_offspring[j]);
                    fitnesses.push_back(t);
                    min = t < min || min == -1 ? t : min;
                }
            }
        }
        population.clear();
        population = new_population;
        k++;
    }
    return min;
}
