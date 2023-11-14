#include <iostream>
#include <vector>
#include <cmath>

typedef struct vehicle
{
    int capacity;
    int time;
    int x_cord;
    int y_cord;
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
} receiver;

int time_to_arrive(vehicle v, receiver r);
void heapify(std::vector<receiver> receivers, std::vector<double> coefficient, int n, int i);
void heapsort(std::vector<receiver> receivers, std::vector<double> coefficient, int n);
int greedy(std::vector<vehicle> vehicles, std::vector<receiver> receivers);

int main(int argc, char *argv[])
{
    std::vector<vehicle> vehicles;
    std::vector<receiver> receivers;
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

int greedy(std::vector<vehicle> vehicles, std::vector<receiver> receivers)
{
    std::vector<double> coefficients;
    std::vector<receiver> receivers_copy;
    std::copy(receivers, receivers_copy, receivers.size());
    int i = 0, j = 0;
    while (i < receivers_copy.size())
    {
        for (int k = 0; k < receivers_copy.size(); k++)
        {
            coefficients.push_back(time_to_arrive(vehicles[j], receivers_copy[k]) / receivers_copy[k].demand);
        }
        heapsort(receivers_copy, coefficients, receivers_copy.size());
        int t = time_to_arrive(vehicles[j], receivers_copy[i]);
        if (vehicles[j].capacity >= receivers_copy[i].demand &&
            vehicles[j].time + t <= receivers_copy[i].time_window_end)
        {
            vehicles[j].capacity -= receivers_copy[i].demand;
            vehicles[j].time += t + receivers_copy[i].service_time;
            vehicles[j].x_cord = receivers_copy[i].x_cord;
            vehicles[j].y_cord = receivers_copy[i].y_cord;
            receivers_copy.erase(receivers_copy.begin() + i);
            coefficients.clear();
            i++;
        }
        else
        {
            j++;
            if (j == vehicles.size())
            {
                return 0;
            }
            continue;
        }
    }
    int result = 0;
    for (auto vehicle : vehicles)
    {
        result += vehicle.time;
    }
    return result;
}