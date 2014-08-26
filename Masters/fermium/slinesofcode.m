%count lines of code
clear all
close all
clc

names=textread('filenames.txt','%s');
linesofcode=0;
for i=1:length(names);
    linesofcode=linesofcode+sloc(names{i});
end
linesofcode