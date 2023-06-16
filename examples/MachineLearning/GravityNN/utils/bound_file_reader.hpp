#pragma once

#include <network/types.hpp>

#include <fstream>
#include <sstream>
#include <vector>
#include <regex>

std::vector<Bound> readBoundsFromFile(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<Bound> bounds;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);

        std::string layer_and_neuron;
        std::getline(iss, layer_and_neuron, ':');
        std::size_t comma_pos = layer_and_neuron.find(',');
        std::string layer_name = layer_and_neuron.substr(0, comma_pos);
        std::string neuron_name = layer_and_neuron.substr(comma_pos + 1);

        std::string bound;
        std::getline(iss, bound);
        bound = bound.substr(2, bound.size() - 3); // Removes " [" at the start and "]" at the end.

        std::regex rgx("\\s*,\\s*");
        std::sregex_token_iterator iter(bound.begin(), bound.end(), rgx, -1);
        std::sregex_token_iterator end;
        std::vector<std::string> splits(iter, end);

        float lower_bound = std::stof(splits[0]);
        float upper_bound = std::stof(splits[1]);

        bounds.push_back(Bound(layer_name, neuron_name, lower_bound, LOWER));
        bounds.push_back(Bound(layer_name, neuron_name, upper_bound, UPPER));
    }

    return bounds;
}
