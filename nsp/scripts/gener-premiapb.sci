P = premia_create();
P.set_asset[str="equity"];
P.set_model[str="BlackScholesndim"];

// les options sont réparties en famille
// premia_get_family (famille, model) renvoie les options de la famille qui
// existe dans le modèle
P.set_option[8, 3];

// premia_get_methods (famille, option, modèle) renvoie la liste des méthodes
// utilisables
P.set_method[1] // MC_LongstaffSchwatrz_ND

p_opt = P.get_option_values[];

p_opt(2)(3) = 2; // met la mautirté à 2
p_opt(3)(3) = 120; // met le strike mautirté à 120

P.set_option_values[p_opt]; // utilise la liste p_opt pour définir les
                            // paramètres de l'option

// les mêmes fonctions existent avec model et method à la place d'option

save('pipo',P); // sauve le pb dans pipo
