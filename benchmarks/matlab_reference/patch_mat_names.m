args = argv(); src = args{1}; crnp = args{2}; dst = args{3};
load(src);                                   % Results_balanced class_with_balanced B CC group
crn = jsondecode(fileread(crnp));
m = Results_balanced.MODEL_r{1};
m.mets = crn.metabolites(:); m.metNames = crn.metabolites(:); m.rxns = crn.reactions(:);
Results_balanced.MODEL_r{1} = m;
save('-v7', dst, 'Results_balanced', 'class_with_balanced', 'B', 'CC', 'group');
printf('%s: fields now %s\n', dst, strjoin(fieldnames(Results_balanced.MODEL_r{1})', ' '));
