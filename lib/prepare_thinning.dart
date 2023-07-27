import 'dart:math';

import 'package:collection/collection.dart';

List<List<List<double>>> prepareThinning(
    {List<Map<String, dynamic>>? thinning,
    List<String> spNames = const ['Fagus sylvatica', 'Pinus sylvestris']}) {
  if (thinning == null) {
    return List.generate(
        1,
        (_) => List.generate(
            5, (_) => List.generate(spNames.length, (_) => double.nan)));
  } else {
    final columnNames = [
      'species',
      'age',
      'stems_n',
      'stem',
      'root',
      'foliage'
    ];
    final thinningData = thinning
        .where((row) => spNames.contains(row['species']))
        .map((row) => Map.fromEntries(
            row.entries.where((entry) => columnNames.contains(entry.key))))
        .toList();
    if (thinningData.isEmpty) {
      return List.generate(
          1,
          (_) => List.generate(
              5, (_) => List.generate(spNames.length, (_) => double.nan)));
    } else {
      thinningData
          .sort((a, b) => (a['species'] as String).compareTo(b['species']));
      thinningData.sort((a, b) => (a['age'] as int).compareTo(b['age'] as int));

      final t_t = Map<String, int>.fromIterable(
          thinningData.map((row) => row['species'] as String),
          key: (species) => species,
          value: (species) =>
              thinningData.where((row) => row['species'] == species).length);
      final n_man = t_t.values.fold(0, max);

      final thinningTable = thinningData.expand((row) {
        final species = spNames.indexOf(row['species'] as String) + 1;
        final thinN = Iterable.generate(t_t[row['species']]!, (i) => i + 1);
        final thinningRows =
            thinN.map((n) => {'species': species, 'thin_n': n});
        return thinningRows;
      }).toList();

      final mergedTable =
          groupBy(thinningTable, (row) => [row['species'], row['thin_n']]);
      final thinningArray = mergedTable.entries
          .map((entry) => (entry.value as List)
              .map((row) => [
                    row['stems_n'] as double,
                    row['stem'] as double,
                    row['root'] as double,
                    row['foliage'] as double
                  ])
              .toList())
          .toList();

      return thinningArray;
    }
  }
}
