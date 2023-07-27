import '3pg.dart';

void run_3PG(
  site,
  species,
  climate, {
  thinning,
  parameters,
  size_dist,
  settings,
  bool check_input = true,
  bool df_out = true,
}) {
  bool thinn_null = thinning == null;
  if (check_input) {
    var input_checked = prepare_input(
        site: site,
        species: species,
        climate: climate,
        thinning: thinning,
        parameters: parameters,
        size_dist: size_dist,
        settings: settings);
  }
}

List prepare_input(
  
{required site,
  required species,
  required climate,
  thinning ,
  parameters ,
  size_dist ,
  Settings? settings ,}
){

  // Site
  site = prepare_site(site = site);

  // Species
  species = prepare_species(species = species);

  // Settings
  Settings set_def = settings?? Settings( 1, 1, 1,  1, 0, 0);

  // Climate
  if( set_def.calculate_d13c == 1 ){
    if( !all( c("co2","d13catm") %in% colnames(climate) ) ){
      stop('Please provide forcing data for co2 and d13catm in climate, if calculate_d13c = 1')
    }
  }

  climate = prepare_climate(climate = climate, from = site$from, to = site$to)

  // Thinning
  thinning = prepare_thinning( thinning = thinning, sp_names = species$species)

  // Parameters
  parameters = prepare_parameters( parameters = parameters, sp_names = species$species)

  // Size distribution
  if( set_def['correct_bias'] == 1 & is.null(size_dist) ){
    stop('Please provide size_dist table or change the setting to size_dist = 0')
  }
  size_dist = prepare_sizeDist( size_dist = size_dist, sp_names = species$species)


  // return the checked output
  return ( site = site, species = species, climate = climate, thinning = thinning, parameters = parameters, size_dist = size_dist, settings = set_def)

}
