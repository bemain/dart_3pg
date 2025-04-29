class Species {
  final String species;
  final String planted;
  final double fertility;
  final double stems_n;
  final double biom_stem;
  final double biom_root;
  final double biom_foliage;

  Species({
    required this.species,
    required this.planted,
    required this.fertility,
    required this.stems_n,
    required this.biom_stem,
    required this.biom_root,
    required this.biom_foliage,
  });
}

Species prepareSpecies(Map<String, dynamic> species) {
  if (species.length != 7) {
    throw Exception(
        'Species table must contain columns: species, planted, fertility, stems_n, biom_stem, biom_root, biom_foliage');
  }

  if (!species.containsKey('species') ||
      !species.containsKey('planted') ||
      !species.containsKey('fertility') ||
      !species.containsKey('stems_n') ||
      !species.containsKey('biom_stem') ||
      !species.containsKey('biom_root') ||
      !species.containsKey('biom_foliage')) {
    throw Exception(
        'Columns names of the species table must correspond to: species, planted, fertility, stems_n, biom_stem, biom_root, biom_foliage');
  }

  // Test for NA
  if (species.containsValue(null)) {
    throw Exception('Species table should not contain null values');
  }

  final speciesName = species['species'].toString();
  final planted = species['planted'].toString();
  final fertility = species['fertility'].toDouble();
  final stems_n = species['stems_n'].toDouble();
  final biom_stem = species['biom_stem'].toDouble();
  final biom_root = species['biom_root'].toDouble();
  final biom_foliage = species['biom_foliage'].toDouble();

  if (fertility < 0 || fertility > 1) {
    throw Exception('Fertility shall be within a range of [0:1]');
  }

  if (stems_n < 0) {
    throw Exception('Stem number shall be greater than 0');
  }

  if (biom_stem < 0) {
    throw Exception('Biomass stem shall be greater than 0');
  }

  if (biom_root < 0) {
    throw Exception('Biomass root shall be greater than 0');
  }

  if (biom_foliage < 0) {
    throw Exception('Biomass foliage shall be greater than 0');
  }

  return Species(
    species: speciesName,
    planted: planted,
    fertility: fertility,
    stems_n: stems_n,
    biom_stem: biom_stem,
    biom_root: biom_root,
    biom_foliage: biom_foliage,
  );
}
