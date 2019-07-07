%Convert date to julian day
%Simple but useful
%By Félix Comtois-Boutet
%Created 23 january 2005
function [julian_day]=date2julien(year,month,day)

%Specify day 0 
number_of_day_before_day_one=datenum(1980,01,00); %Number of day before day one on IBM-PC
absolute_number_of_day=datenum(year,month,day); %Default day one for datenum is january 1st of year 0
julian_day=absolute_number_of_day - number_of_day_before_day_one; %Number of day since day one