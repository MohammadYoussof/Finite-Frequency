function rayxyz = go_project_xy( par, raylat, raylon )

z1 = par.vel.z1;
z2 = par.vel.z2;

[rayx,rayy] = project_xy( par, raylat, raylon );
rayz = [z1(1); z2];
rayxyz = [rayx rayy rayz];
