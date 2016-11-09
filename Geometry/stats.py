def get_stats(redcoe_obj, percent=True):
    stats_dict = {'strong': {}, 'weak': {}}
    event_count = len(redcoe_obj._dict)
    for i in range(1, redcoe_obj._range + 1):
        strong_sum = 0
        weak_sum = 0
        for event in redcoe_obj._dict:  # event is a node in a triangle
            strong_sum += redcoe_obj._dict[event]['strong'][i]
            weak_sum += redcoe_obj._dict[event]['weak'][i]
        if percent:
            try:
                stats_dict['strong'][i] = strong_sum / event_count
            except ZeroDivisionError:  # dict is empty
                raise Exception("Graph has not been evaluated yet")
            stats_dict['weak'][i] = weak_sum / event_count
        else:
            stats_dict['strong'][i] = strong_sum
            stats_dict['weak'][i] = weak_sum
    return stats_dict
