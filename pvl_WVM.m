function [smooth_irradiance,other_outputs]=pvl_WVM(irr_sensor,plantinfo,cloud_speed,otherinputs)
%PVL_WVM Wavelet Variability Model
%
%Syntax:
% smooth_irradiance = pvl_wvm(irr_sensor,plantinfo,cloud_speed) computes the
%   spatially-smoothed irradiance using the wavelet variability model. 
%
%irr_sensor is a struct with variables:
%   irr_sensor.irr: the irradiance measurement
%   irr_sensor.time: the time stamps (Matlab time vector) for irr_sensor.irr
%   irr_sensor.Lat: latitude of the sensor
%   irr_sensor.Lon: longitude of the sensor
%   irr_sensor.alt: altitude of the sensor
%   irr_sensor.tilt: tilt angle of the sensor, 0 = flat (e.g., GHI)
%   irr_sensor.azimuth: azimuth angle of the sensor, 180 = due south
%   (optional) irr.sensor.clear_sky_irradiance:  manually enter the clear-sky irradiance (e.g., for an irradiance sensor on a tracking system) 
%   irr_sensor.UTCoffset: UTC offset
%
%plantinfo is a struct describing the plant to simulate with variables:
%   plantinfo.tilt: tilt angle of plant modules 
%   plantinfo.azimuth: azimuth angle of plant modules
%   (optional) plantinfo.clear_sky_irrPOA: manually enter the clear-sky irradiance in the module POA (e.g., for tracking systems)
%   plantinfo.type: 'square','polygon',' or 'discrete'
%       'square' square PV plant with specified number of MWs and PV density
%       'polygon' custom PV plant shape (define vetiticies in lat/lon)
%       'discrete' simulate only certain points (e.g., to replicate output of multiple point sensors)
%   plantinfo.MW: = MW of PV installed (not necessary for 'discrete' type)
%   plantinfo.PVdensity: = W PV installed per m2 in plant area (e.g., 41 W/m2 is 1MW per 6 acres) (not necessary for 'discrete' type)
%   plantinfo.Lat: (only needed for type 'polygon' or 'discrete') latitude of polygon verticies or discrete points
%   plantinfo.Lon: (only needed for type 'polygon' or 'discrete') longitude of polygon verticies or discrete points
%
%cloud_speed is a single value of the daily cloud speed
%
% BSD License:
% Copyright (c) 2012, The Regents of the University of California
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 
% Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 
% Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 
% Neither the name The Regents of the University of California, the names of its campuses nor any abbreviation thereof, nor the names of the contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
p1=which('pvl_WVM');
path1=p1(1:end-length('pvl_WVM.m'));
addpath([path1 'supporting_programs\']);
addpath([path1 'supporting_programs\nansuite']);

% addpath('.\supporting_programs\');
% addpath('.\supporting_programs\nansuite\');

try
    temp=plantinfo.Lon;
catch
    plantinfo.Lon=irr_sensor.Lon;
    plantinfo.Lat=irr_sensor.Lat;
end

%% compute the clear-sky index
if isfield(irr_sensor,'clear_sky_irradiance')==1
    irr_sensor.tilt=NaN;
    irr_sensor.azimuth=NaN;
    [clear_sky_index]=pvl_WVM_compute_clear_sky_index(irr_sensor.time,irr_sensor.irr,irr_sensor.Lat,irr_sensor.Lon,irr_sensor.alt,irr_sensor.UTCoffset,irr_sensor.tilt,irr_sensor.azimuth,irr_sensor.clear_sky_irradiance);
else
    [clear_sky_index]=pvl_WVM_compute_clear_sky_index(irr_sensor.time,irr_sensor.irr,irr_sensor.Lat,irr_sensor.Lon,irr_sensor.alt,irr_sensor.UTCoffset,irr_sensor.tilt,irr_sensor.azimuth);
end
%% compute the wavelet modes
[wavelet,timeout,tmscales]=pvl_WVM_compute_wavelet(irr_sensor.time,clear_sky_index);

%% compute variablity reduction
try
    dist=otherinputs.dist;
catch
if strcmp(plantinfo.type,'discrete')
    [dist]=pvl_WVM_compute_distances(plantinfo.Lon,plantinfo.Lat,plantinfo.type);
end
if strcmp(plantinfo.type,'square')
    [dist]=pvl_WVM_compute_distances(plantinfo.Lon,plantinfo.Lat,plantinfo.type,plantinfo.MW,plantinfo.PVdensity);
end
if strcmp(plantinfo.type,'polygon')
    [dist]=pvl_WVM_compute_distances(plantinfo.Lon,plantinfo.Lat,plantinfo.type,plantinfo.MW);
end
end

VR=pvl_WVM_compute_VR(dist,tmscales,cloud_speed);

%% smooth wavelets by VR 
for i=1:length(tmscales)
    if i<length(tmscales) %special treatment for last timescale to ensure wavelet modes can be recombined to create simulated wavelet
        wavelet_smooth(i,:)=wavelet(i,:)./sqrt(VR(i));
    else
        wavelet_smooth(i,:)=(wavelet(i,:));
    end
end

%% sum wavelets to creat smoothed clear-sky index
 
[C,ia,ib]=intersect(irr_sensor.time,timeout);

clear_sky_index_smooth=zeros(size(clear_sky_index));
clear_sky_index_smooth(ia)=nansum(wavelet_smooth);

%% compute clear-sky irradiance for plant
irradiance_in=ones(size(irr_sensor.time)).*NaN; %temporary varaible needed in running pvl_WVM_compute_clear_sky_index
if isfield(plantinfo,'clear_sky_irrPOA')==1
    clear_sky_irradiance_smooth=plantinfo.clear_sky_irrPOA;
else
    [clear_sky_indexPOA,clear_sky_irradiance_smooth]=pvl_WVM_compute_clear_sky_index(irr_sensor.time,irradiance_in,mean(plantinfo.Lat),mean(plantinfo.Lon),irr_sensor.alt,irr_sensor.UTCoffset,plantinfo.tilt,plantinfo.azimuth);
end
 
%% combine clear-sky index with clear-sky irradiance to generate smoothed output
try
smooth_irradiance=clear_sky_irradiance_smooth.*clear_sky_index_smooth;
catch
    smooth_irradiance=clear_sky_irradiance_smooth.*clear_sky_index_smooth';
end
%% produce other outputs
other_outputs.clear_sky_index_smooth=clear_sky_index_smooth;
other_outputs.VR=VR;
other_outputs.clear_sky_index=clear_sky_index;
other_outputs.wavelet=wavelet;
other_outputs.wavelet_smooth=wavelet_smooth;
other_outputs.time_temp=timeout;
other_outputs.tmscales=tmscales;
other_outputs.dist=dist;
other_outputs.clear_sky_irr_POA=clear_sky_irradiance_smooth;
other_outputs.time=irr_sensor.time;

