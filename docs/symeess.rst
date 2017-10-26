.. highlight:: rst

Symeess
=======
This is the main module

.. automodule:: symeess
   :members:
   :show-inheritance:

.. py:class:: symeess.Symeess

   .. py:method:: write_shape_measure_2file(shape_label, central_atom=None, output_name='../examples/symeess_shape')
      Method that prints to file shape's measure

      :param str shape_label: reference polyhedra label which user will compare with his polyhedra.
                              Reference labels can be found in [#f1]_
      :param int central_atom: position of the central atom in molecule if exist
      :param str output_name: custom name without extension
      :return: shape's measure in the output_name.tab file
   .. py:method:: write_shape_structure_2file(shape_label, central_atom=None, output_name='../examples/symeess_shape')
      Method that prints to file shape's structure

      :param str shape_label: reference polyhedra label which user will compare with his polyhedra.
                              Reference labels can be found in [#f1]_
      :param int central_atom: position of the central atom in molecule if exist
      :param str output_name: custom name without extension
      :return: shape's structure in the output_name.out file


.. [#f1] :doc:`shape_references`