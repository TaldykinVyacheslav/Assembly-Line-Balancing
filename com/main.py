__author__ = 'zegoline'

from functools import reduce
from random import random
from random import choice


def find_path(graph, start, end, path=None):
    if not path:
        path = []
    path = path + [start]
    if start == end:
        return path
    for node in graph[start - 1]:
        if node not in path:
            new_path = find_path(graph, node, end, path)
            if new_path:
                return new_path
    return None


def find_candidates(graph, operations_left, current_time, cycle_time, times_list):
    schedulable_operations = []

    # Find schedulable operations
    for operation_j in operations_left:  # iterate through
        is_operation_candidate = True
        for operation_i in [operation for operation in operations_left if operation != operation_j]:
            if find_path(graph, operation_i, operation_j) is not None:
                is_operation_candidate = False
        if is_operation_candidate:
            schedulable_operations.append(operation_j)

    # Find candidate operations
    candidate_operations = []
    for operation in schedulable_operations:
        if times_list[operation - 1] <= (cycle_time - current_time):
            candidate_operations.append(operation)
    return candidate_operations


# find operations that can be moved to left station
def operations_to_left(operations_left, operations_right, times_list, graph, cycle_time):
    movable_operations = []
    left_station_time = reduce(lambda prev_sum, op: prev_sum + times_list[op - 1], operations_left, 0)

    for operation in operations_right:
        if left_station_time + times_list[operation - 1] > cycle_time:
            continue
        can_move = True

        # Check if candidate operation has no predecessors in right station
        for operation_predecessor in [op_pre for op_pre in operations_right if op_pre != operation]:
            if find_path(graph, operation_predecessor, operation):
                can_move = False
        if can_move:
            movable_operations.append(operation)
    return movable_operations


# find operations that can be moved to right station
def operations_to_right(operations_left, operations_right, times_list, graph, cycle_time):
    movable_operations = []
    right_station_time = reduce(lambda prev_sum, op: prev_sum + times_list[op - 1], operations_right, 0)

    for operation in operations_left:
        if right_station_time + times_list[operation - 1] > cycle_time:
            continue
        can_move = True

        # Check if candidate operation has no successors in left station
        for operation_successor in [op_suc for op_suc in operations_left if op_suc != operation]:
            if find_path(graph, operation, operation_successor):
                can_move = False
        if can_move:
            movable_operations.append(operation)
    return movable_operations


def heuristic_procedure(operations_num, adjacency_matrix, cycle_time, times_list, priority_list):
    operations_left = [op for op in range(1, operations_num + 1)]
    operations_processed = []
    stations_operations = [[]]
    stations_num = 1
    current_time = 0

    while operations_left:
        candidate_operations = find_candidates(adjacency_matrix, operations_left, current_time, cycle_time, times_list)
        # if candidate operations list is empty then add new station
        if not candidate_operations:
            stations_num += 1
            stations_operations.append([])
            current_time = 0
            continue

        # sort by execution time [priority] in descending order
        candidate_operations = sorted(candidate_operations, key=lambda operation: priority_list[operation - 1])
        selected_operation = candidate_operations[0]
        operations_left.remove(selected_operation)
        operations_processed.append(selected_operation)
        current_time += times_list[selected_operation - 1]
        stations_operations[-1].append(selected_operation)
    return stations_num, stations_operations


def local_search_procedure(stations_num, stations_operations, times_list, adjacency_matrix, cycle_time):
    for station in range(1, stations_num):
        while True:
            dirty_flag = False  # flag that indicates that any changes were made
            left_movable_operations = operations_to_left(stations_operations[station - 1], stations_operations[station],
                                                         times_list, adjacency_matrix, cycle_time)
            left_movable_operations = sorted(left_movable_operations, key=lambda op: times_list[op - 1])

            # move operations to left station while left station time does not exceed cycle time
            left_station_time = reduce(lambda prev_sum, op: prev_sum + times_list[op - 1],
                                       stations_operations[station - 1], 0)
            right_station_time = reduce(lambda prev_sum, op: prev_sum + times_list[op - 1],
                                        stations_operations[station], 0)
            for operation in left_movable_operations:
                if left_station_time + times_list[operation - 1] <= cycle_time:
                    stations_operations[station - 1].append(operation)
                    stations_operations[station].remove(operation)
                    left_station_time += times_list[operation - 1]
                    dirty_flag = True
                else:
                    break
            left_movable_operations = operations_to_left(stations_operations[station - 1], stations_operations[station],
                                                         times_list, adjacency_matrix, cycle_time)
            right_movable_operations = operations_to_right(stations_operations[station - 1],
                                                           stations_operations[station], times_list, adjacency_matrix,
                                                           cycle_time)
            exchangeable_operations = [[left_op, right_op] for left_op in left_movable_operations for right_op in
                                       right_movable_operations]
            exchangeable_operations = filter(
                lambda op_pair: (times_list[op_pair[1] - 1] - times_list[op_pair[0] - 1]) > 0, exchangeable_operations)
            exchangeable_operations = sorted(exchangeable_operations,
                                             key=lambda op_pair: times_list[op_pair[1] - 1] - times_list[
                                                 op_pair[0] - 1])
            for operations_pair in exchangeable_operations:
                if left_station_time - times_list[operations_pair[0] - 1] + times_list[
                            operations_pair[1] - 1] <= cycle_time:
                    # move operations_pair[1] to left station
                    stations_operations[station - 1].append(operations_pair[1])
                    stations_operations[station].remove(operations_pair[1])
                    # move operations_pair[0] to right station
                    stations_operations[station].append(operations_pair[0])
                    stations_operations[station - 1].remove(operations_pair[0])
                    # recalculate left and right stations times
                    left_station_time += (times_list[operations_pair[1] - 1] - times_list[operations_pair[0] - 1])
                    right_station_time += (times_list[operations_pair[0] - 1] - times_list[operations_pair[1] - 1])
                    dirty_flag = True
                else:
                    break

            if not dirty_flag:
                break

    # try to reduce number of work stations
    while not stations_operations[-1]:
        del stations_operations[-1]
        stations_num -= 1
    return stations_num, stations_operations


def genetic_algorithm_procedure(population_num, generations_num, operations_num, adjacency_matrix, cycle_time,
                                times_list, priority_list, stations_num, stations_operations, elitist_part,
                                coin_probability, mutation_probability):
    population = []
    fitness_func_values = {}

    # create first generation
    for i in range(1, population_num):
        population.append([random() for _ in range(operations_num)])
        # calculate fitness function values for each chromosome in population
        for chromosome in population:
            # Heuristic priority-based procedure
            (chromosome_stations_num, chromosome_stations_ops) = heuristic_procedure(operations_num, adjacency_matrix,
                                                                                     cycle_time, times_list,
                                                                                     priority_list)
            # The local search procedure
            (chromosome_stations_num, chromosome_stations_ops) = local_search_procedure(stations_num,
                                                                                        chromosome_stations_ops,
                                                                                        times_list, adjacency_matrix,
                                                                                        cycle_time)
            fitness_func_values[''.join([('%.5f' % x) for x in chromosome])] = {'stations_num': chromosome_stations_num,
                                                                                'stations_ops': chromosome_stations_ops}

    # iterations for each generation
    for generation_no in range(generations_num):
        # reproduction stage
        new_population = list(population[0: int(elitist_part * len(population))])
        # crossovers stage
        for i in range(int(population_num * (1 - elitist_part - mutation_probability))):
            parent1, parent2 = choice(population), choice(population)
            new_gene = []
            for gene1, gene2 in zip(parent1, parent2):
                if random() < coin_probability:
                    new_gene.append(gene1)
                else:
                    new_gene.append(gene2)
            new_population.append(new_gene)
        # mutation stage
        for i in range(int(population_num * mutation_probability)):
            new_population.append([random() for _ in range(operations_num)])
        population = new_population

        # calculate fitness function values for each chromosome in population
        for chromosome in population:
            # Heuristic priority-based procedure
            (chromosome_stations_num, chromosome_stations_ops) = heuristic_procedure(operations_num, adjacency_matrix,
                                                                                     cycle_time, times_list,
                                                                                     priority_list)
            # The local search procedure
            (chromosome_stations_num, chromosome_stations_ops) = local_search_procedure(stations_num,
                                                                                        chromosome_stations_ops,
                                                                                        times_list, adjacency_matrix,
                                                                                        cycle_time)
            fitness_func_values[''.join([('%.5f' % x) for x in chromosome])] = {'stations_num': chromosome_stations_num,
                                                                                'stations_ops': chromosome_stations_ops}
        # sort population by fitness function values
        population.sort(
            key=lambda chromosome: fitness_func_values[''.join([('%.5f' % x) for x in chromosome])]['stations_num'])

    # get the best individual from the last population and compare its value with the found value
    best_value = fitness_func_values[''.join([('%.5f' % x) for x in population[0]])]
    if best_value['stations_num'] < stations_num:
        stations_num = best_value['stations_num']
        stations_operations = best_value['stations_ops']
    return stations_num, stations_operations


def main():
    # Read input file
    # Input file structure:
    # 1. Cycle time
    # 2. Operations number
    # 3. List of execution times for each operation
    # 4. Matrix of dependencies of operations
    # 5. Population size for genetic algorithm
    # 6. Number of generations
    # 7. Part of population that will be mutated
    # 8. Part of population that will be copied to next generation in elitist strategy
    # 9. Coin probability used to choose specific gene from one of parent chromosomes
    with open('../resources/input.txt') as f:
        cycle_time = int(f.readline())
        operations_num = int(f.readline())
        times_list = [(int(x)) for x in f.readline().split()]
        priority_list = [(int(x)) for x in times_list]
        priority_list = [x / max(priority_list) for x in priority_list]
        adjacency_matrix = [[(int(x)) for x in f.readline().split()] for _ in range(operations_num)]
        # genetic algorithm parameters
        population_num = int(f.readline())
        generations_num = int(f.readline())
        mutation_probability = float(f.readline())
        elitist_part = float(f.readline())
        coin_probability = float(f.readline())

    # Heuristic priority-based procedure
    (stations_num, stations_operations) = heuristic_procedure(operations_num, adjacency_matrix, cycle_time, times_list,
                                                              priority_list)

    # The local search procedure
    (stations_num, stations_operations) = local_search_procedure(stations_num, stations_operations, times_list,
                                                                 adjacency_matrix, cycle_time)

    # Genetic algorithm and chromosome representation
    (stations_num, stations_operations) = genetic_algorithm_procedure(population_num, generations_num, operations_num,
                                                                      adjacency_matrix, cycle_time, times_list,
                                                                      priority_list, stations_num, stations_operations,
                                                                      elitist_part, coin_probability,
                                                                      mutation_probability)

    # Write results
    # Output file structure:
    # 1. Final minimal number of stations
    # 2. Number of operations
    # 3. Lists of operations per each station
    with open("../resources/output.txt", "w") as f:
        f.write(str(stations_num) + '\n')
        f.write(str(operations_num) + '\n')
        for station_operations in stations_operations:
            f.write(' '.join([str(x) for x in station_operations]) + '\n')
    print('Station number is %d' % stations_num)


if __name__ == "__main__":
    main()
