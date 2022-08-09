"""
Miscellaneous tool functions for use elsewhere

"""

from numbers import Number
import numpy as np
import matplotlib.pyplot as plt

def is_numeric(x, positive=False, not_negative=False, unit=False):
    """
    Check if input x is either numeric, or is a list of numeric
    """
    if isinstance(x, Number):
        if positive:
            # x must be greater than zero
            if x > 0:
                return(True)
            else:
                return(False)

        elif not_negative:
            # x must be greater than or equal to zero
            if x >= 0:
                return(True)
            else:
                return(False)

        elif unit:
            # x must be between 0 and 1
            if x>=0 and x<=1:
                return(True)
            else:
                return(False)

        else:
            return(True)

    elif type(x) == tuple or type(x) == list:
        all_items_numeric = all(is_numeric(item, positive, not_negative) for item in x)
        return(all_items_numeric)
    
    else:
        return(False)

def plot_stacked_bar(data, series_labels, category_labels=None, 
                     show_values=False, value_format="{}", y_label=None, 
                     colors=None, grid=True, reverse=False):
    """Plots a stacked bar chart with the data and labels provided.

    Keyword arguments:
    data            -- 2-dimensional numpy array or nested list
                       containing data for each series in rows
    series_labels   -- list of series labels (these appear in
                       the legend)
    category_labels -- list of category labels (these appear
                       on the x-axis)
    show_values     -- If True then numeric value labels will 
                       be shown on each bar
    value_format    -- Format string for numeric value labels
                       (default is "{}")
    y_label         -- Label for y-axis (str)
    colors          -- List of color labels
    grid            -- If True display grid
    reverse         -- If True reverse the order that the
                       series are displayed (left-to-right
                       or right-to-left)
    """

    ny = len(data[0])
    ind = list(range(ny))

    axes = []
    cum_size = np.zeros(ny)

    data = np.array(data)

    if reverse:
        data = np.flip(data, axis=1)
        category_labels = reversed(category_labels)

    for i, row_data in enumerate(data):
        color = colors[i] if colors is not None else None
        axes.append(plt.bar(ind, row_data, bottom=cum_size, 
                            label=series_labels[i], color=color))
        cum_size += row_data

    if category_labels:
        plt.xticks(ind, category_labels)

    if y_label:
        plt.ylabel(y_label)

    plt.legend()

    if grid:
        plt.grid()

    if show_values:
        for axis in axes:
            for bar in axis:
                w, h = bar.get_width(), bar.get_height()
                plt.text(bar.get_x() + w/2, bar.get_y() + h/2, 
                         value_format.format(h), ha="center", 
                         va="center")

def plot_power_balance(heat_reports, labels):
    """"
    
        heat_reports (list of dicts):    list of standardised reports
    """
    
    gain_keys = []
    gain_values = [ [] for _ in range(len(heat_reports)) ]
    loss_keys = []
    loss_values = [ [] for _ in range(len(heat_reports)) ]

    for key in heat_reports[0].keys():
        if "gain" in key:
            gain_keys.append(key)
            for index, report in enumerate(heat_reports):
                gain_values[index].append(report[key])
        
        elif "loss" in key:
            loss_keys.append(key)
            for index, report in enumerate(heat_reports):
                loss_values[index].append(report[key])
    
    
    for index in range(len(heat_reports)):
        gain_values.insert((2*index)+1, 
            np.zeros(len(gain_values[0])))
        loss_values.insert((2*index), 
            np.zeros(len(loss_values[0])))

    cum_gain = np.zeros(len(heat_reports)*2)
    gain_values = np.asarray(gain_values)
    loss_values = np.asarray(loss_values)
    all_keys = gain_keys + loss_keys

    #labels = ["a", "b"]

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    shifted=False
    
    for value_set in range(len(all_keys)):
        if all_keys[value_set] in gain_keys:
            plt.bar(x-width/2, 
            gain_values[::2,value_set], 
            width, 
            label=all_keys[value_set], 
            bottom=cum_gain[::2])


            for x_item in range(len(x)):
                if gain_values[x_item*2,value_set] != 0:
                    plt.text(x[x_item]-width, 
                        gain_values[x_item*2,value_set]/2 + cum_gain[x_item*2], 
                        all_keys[value_set], 
                        va="center")

            cum_gain = np.add(cum_gain, gain_values[:,value_set])

            

        elif all_keys[value_set] in loss_keys:

            offset_value_set = value_set - len(gain_keys)

            plt.bar(x+width/2, 
            loss_values[1::2,offset_value_set], 
            width, 
            label=all_keys[value_set], 
            bottom=np.add(cum_gain[1::2], cum_gain[::2]))

            for x_item in range(len(x)):
                if loss_values[1 + x_item*2,offset_value_set] != 0:
                    plt.text(x[x_item], 
                        loss_values[1 + x_item*2,offset_value_set]/2 + np.add(cum_gain[1 + x_item*2], cum_gain[x_item*2]), 
                        all_keys[value_set], 
                        va="center")
            
            cum_gain = np.add(cum_gain, loss_values[:,offset_value_set])

        
    plt.xticks(x, labels) 

    plt.ylabel("Heat loss into habitat (W) - positive=loss")   
    
    plt.legend(bbox_to_anchor = (1.05, 0.6))
    
    
    