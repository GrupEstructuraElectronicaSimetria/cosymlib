from symeess.shape import maps


def write_minimum_distortion_path_shape_2file(shape_label1, shape_label2, num_points=20, show=False, output_name=None):
    path = get_shape_map(shape_label1, shape_label2, num_points)
    if show:
        import matplotlib.pyplot as plt
        plt.plot(path[0], path[1], 'k', linewidth=2.0)
        plt.xlabel(shape_label1)
        plt.ylabel(shape_label2)
        plt.show()
    else:
        from symeess.file_io import shape2file
        shape2file.write_shape_map(shape_label1, shape_label2, path, output_name)


def get_shape_map(shape_label1, shape_label2, num_points):
    return maps.get_shape_map(shape_label1, shape_label2, num_points)
