% Implementation of a published MICE model by Morello et al. 2014, MEPS
% Karlo Hock, University of Queensland, v1 2014; v2 2019

% any number of years
years = 18;

% runs an example scenario of the COTS-coral population dynamics model published by Morello et al 2014, MEPS;
% See COTS_MICE_popdyn file for more details;
% cots_numbers returns all three size categories from the model
[cots_numbers] = COTS_MICE_popdyn(years);

% plot the adult population over time
plot(0:18, cots_numbers(:,3));
title('Adult COTS abundance', 'FontSize', 11);
xlabel('Years');
ylabel('COTS per tow');
