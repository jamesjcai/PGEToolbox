function [RR] = hgmutrate(seq)
% Hwang, D. G. & Green, P. (2004) PNAS 101, 13994–14001.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[n, m] = size(seq);
if m < 3,
    error('Too short')
end
[Rz] = i_getRz;

RR = [];

for ss = 1:n
    R = zeros(4, m);
    for k = 1:m - 2
        a = seq(ss, k);
        b = seq(ss, k+1);
        c = seq(ss, k+2);
        if ~any(seq(ss, k:k+2) > 4)
            z = dinuclise16([a, c]);
            for x = 1:4
                R(x, k+1) = Rz(z, dinuclise16([b, x]));
            end
        end
    end
    RR = [RR; R];
end


    function [Rz] = i_getRz

        d = [1, 5, 0.62769, 0.027638; 1, 7, 0.45859, 0.022966; 1, 8, 2.6823, 0.058583; ...
            2, 5, 0.83118, 0.040471; 2, 7, 0.6428, 0.033228; 2, 8, 2.8617, 0.075674; ...
            3, 5, 1.6006, 0.16721; 3, 7, 1.5488, 0.16528; 3, 8, 21.456, 0.68138; ...
            4, 5, 0.58017, 0.028806; 4, 7, 0.59816, 0.02538; 4, 8, 2.6724, 0.060996; ...
            1, 13, 0.58196, 0.024333; 1, 14, 2.1416, 0.047116; 1, 15, 0.52335, 0.022251; ...
            2, 13, 0.45699, 0.027302; 2, 14, 1.6567, 0.050393; 2, 15, 0.45096, 0.025288; ...
            3, 13, 0.53882, 0.026186; 3, 14, 1.9493, 0.056762; 3, 15, 0.61176, 0.026369; ...
            4, 13, 0.41126, 0.017032; 4, 14, 1.6981, 0.035344; 4, 15, 0.43822, 0.016898; ...
            5, 5, 0.4542, 0.026465; 5, 7, 0.2539, 0.019908; 5, 8, 2.1802, 0.051212; ...
            6, 5, 0.59279, 0.035225; 6, 7, 0.41379, 0.029179; 6, 8, 3.1005, 0.072224; ...
            7, 5, 1.8446, 0.18715; 7, 7, 1.6297, 0.18157; 7, 8, 17.357, 0.59105; ...
            8, 5, 0.42671, 0.028342; 8, 7, 0.41135, 0.027173; 8, 8, 2.6498, 0.057512; ...
            5, 13, 0.46528, 0.030403; 5, 14, 1.6726, 0.05398; 5, 15, 0.50977, 0.029703; ...
            6, 13, 0.40515, 0.025687; 6, 14, 1.1593, 0.039525; 6, 15, 0.34929, 0.023798; ...
            7, 13, 0.42327, 0.025102; 7, 14, 1.173, 0.038948; 7, 15, 0.4618, 0.026361; ...
            8, 13, 0.40686, 0.024463; 8, 14, 1.126, 0.039039; 8, 15, 0.46515, 0.02342; ...
            9, 5, 0.82222, 0.034431; 9, 7, 0.36573, 0.024599; 9, 8, 2.1361, 0.059787; ...
            10, 5, 0.69526, 0.038472; 10, 7, 0.54949, 0.033149; 10, 8, 2.6337, 0.075095; ...
            11, 5, 2.0471, 0.2024; 11, 7, 1.0873, 0.15253; 11, 8, 12.803, 0.44905; ...
            12, 5, 0.567, 0.029239; 12, 7, 0.40006, 0.025687; 12, 8, 2.3241, 0.060712; ...
            9, 13, 0.62273, 0.035996; 9, 14, 2.0633, 0.0628; 9, 15, 0.59245, 0.032905; ...
            10, 13, 0.38124, 0.029113; 10, 14, 1.6557, 0.061536; 10, 15, 0.42416, 0.032519; ...
            11, 13, 0.5239, 0.032256; 11, 14, 1.1339, 0.045077; 11, 15, 0.70483, 0.032827; ...
            12, 13, 0.3634, 0.023244; 12, 14, 1.7942, 0.059134; 12, 15, 0.47632, 0.025693; ...
            13, 5, 0.49111, 0.023168; 13, 7, 0.37711, 0.019188; 13, 8, 1.8721, 0.042589; ...
            14, 5, 0.45557, 0.02522; 14, 7, 0.48697, 0.02558; 14, 8, 2.6698, 0.058641; ...
            15, 5, 1.9155, 0.19839; 15, 7, 2.2981, 0.20717; 15, 8, 17.398, 0.58224; ...
            16, 5, 0.60825, 0.023714; 16, 7, 0.5799, 0.023643; 16, 8, 1.8472, 0.040011; ...
            13, 13, 0.5641, 0.020881; 13, 14, 1.1798, 0.032071; 13, 15, 0.58885, 0.024382; ...
            14, 13, 0.23936, 0.016064; 14, 14, 1.1621, 0.037588; 14, 15, 0.36264, 0.019927; ...
            15, 13, 0.28258, 0.018383; 15, 14, 1.0884, 0.042083; 15, 15, 0.67395, 0.026783; ...
            16, 13, 0.29689, 0.013341; 16, 14, 0.86382, 0.025514; 16, 15, 0.48529, 0.017936; ...
            16, 12, 0.62769, 0.027638; 16, 10, 0.45859, 0.022966; 16, 9, 2.6823, 0.058583; ...
            12, 12, 0.83118, 0.040471; 12, 10, 0.6428, 0.033228; 12, 9, 2.8617, 0.075674; ...
            8, 12, 1.6006, 0.16721; 8, 10, 1.5488, 0.16528; 8, 9, 21.456, 0.68138; ...
            4, 12, 0.58017, 0.028806; 4, 10, 0.59816, 0.02538; 4, 9, 2.6724, 0.060996; ...
            16, 4, 0.58196, 0.024333; 16, 3, 2.1416, 0.047116; 16, 2, 0.52335, 0.022251; ...
            12, 4, 0.45699, 0.027302; 12, 3, 1.6567, 0.050393; 12, 2, 0.45096, 0.025288; ...
            8, 4, 0.53882, 0.026186; 8, 3, 1.9493, 0.056762; 8, 2, 0.61176, 0.026369; ...
            4, 4, 0.41126, 0.017032; 4, 3, 1.6981, 0.035344; 4, 2, 0.43822, 0.016898; ...
            15, 12, 0.4542, 0.026465; 15, 10, 0.2539, 0.019908; 15, 9, 2.1802, 0.051212; ...
            11, 12, 0.59279, 0.035225; 11, 10, 0.41379, 0.029179; 11, 9, 3.1005, 0.072224; ...
            7, 12, 1.8446, 0.18715; 7, 10, 1.6297, 0.18157; 7, 9, 17.357, 0.59105; ...
            3, 12, 0.42671, 0.028342; 3, 10, 0.41135, 0.027173; 3, 9, 2.6498, 0.057512; ...
            15, 4, 0.46528, 0.030403; 15, 3, 1.6726, 0.05398; 15, 2, 0.50977, 0.029703; ...
            11, 4, 0.40515, 0.025687; 11, 3, 1.1593, 0.039525; 11, 2, 0.34929, 0.023798; ...
            7, 4, 0.42327, 0.025102; 7, 3, 1.173, 0.038948; 7, 2, 0.4618, 0.026361; ...
            3, 4, 0.40686, 0.024463; 3, 3, 1.126, 0.039039; 3, 2, 0.46515, 0.02342; ...
            14, 12, 0.82222, 0.034431; 14, 10, 0.36573, 0.024599; 14, 9, 2.1361, 0.059787; ...
            10, 12, 0.69526, 0.038472; 10, 10, 0.54949, 0.033149; 10, 9, 2.6337, 0.075095; ...
            6, 12, 2.0471, 0.2024; 6, 10, 1.0873, 0.15253; 6, 9, 12.803, 0.44905; ...
            2, 12, 0.567, 0.029239; 2, 10, 0.40006, 0.025687; 2, 9, 2.3241, 0.060712; ...
            14, 4, 0.62273, 0.035996; 14, 3, 2.0633, 0.0628; 14, 2, 0.59245, 0.032905; ...
            10, 4, 0.38124, 0.029113; 10, 3, 1.6557, 0.061536; 10, 2, 0.42416, 0.032519; ...
            6, 4, 0.5239, 0.032256; 6, 3, 1.1339, 0.045077; 6, 2, 0.70483, 0.032827; ...
            2, 4, 0.3634, 0.023244; 2, 3, 1.7942, 0.059134; 2, 2, 0.47632, 0.025693; ...
            13, 12, 0.49111, 0.023168; 13, 10, 0.37711, 0.019188; 13, 9, 1.8721, 0.042589; ...
            9, 12, 0.45557, 0.02522; 9, 10, 0.48697, 0.02558; 9, 9, 2.6698, 0.058641; ...
            5, 12, 1.9155, 0.19839; 5, 10, 2.2981, 0.20717; 5, 9, 17.398, 0.58224; ...
            1, 12, 0.60825, 0.023714; 1, 10, 0.5799, 0.023643; 1, 9, 1.8472, 0.040011; ...
            13, 4, 0.5641, 0.020881; 13, 3, 1.1798, 0.032071; 13, 2, 0.58885, 0.024382; ...
            9, 4, 0.23936, 0.016064; 9, 3, 1.1621, 0.037588; 9, 2, 0.36264, 0.019927; ...
            5, 4, 0.28258, 0.018383; 5, 3, 1.0884, 0.042083; 5, 2, 0.67395, 0.026783; ...
            1, 4, 0.29689, 0.013341; 1, 3, 0.86382, 0.025514; 1, 2, 0.48529, 0.017936];

        Rz = zeros(16);
        for k = 1:size(d, 1)
            Rz(d(k, 1), d(k, 2)) = d(k, 3);
        end

        DINU = {'AA', 'AC', 'AG', 'AT', ...
            'CA', 'CC', 'CG', 'CT', ...
            'GA', 'GC', 'GG', 'GT', ...
            'TA', 'TC', 'TG', 'TT'};
        % matrixcircle(Rz,DINU);
