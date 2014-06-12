classdef SuperModel < handle
    % SUPERMODEL base model which all other models extend
    
    properties
     
        stimulusLoader
        
        stimulusSize
        
    end
    
    properties (Constant)
        
        DT = 0.001;
        
    end
    
    methods

        function SM = SuperModel()
        end
        
        function OUT = stimulate(SM, stimulus)
        % STIMULATE stimulates the model using a stimulus
        % 
        % PARAMS
        %   (optional: if none provided, uses stimulus_loader) stimulus 
        %
        % OUTPUT
        %   OUT.firing_rate
        %   OUT.spike        
        end
        
        function [] = setStimulusLoader(SM, loader)
        % SETSTIMULUSLOADER saves a stimulus loader and updates the kernel
        % to match the stimulus size (and other relevant parameters)
            SM.stimulusLoader = loader;
            stimulus = SM.stimulusLoader.randomStimulus();
            SM.stimulusSize = size(stimulus);
        end
        
        function [] = resetCache(SM)
        end
        
        function [] = fit(SM)
        end
        
        function [] = plotKernel(SM)
        end
    end
    
    methods(Static)
        
        function [spike, noisyRate] = rate2spike(rate, dt)
        % RATE2SPIKE adds noise to firing RATE, and then generates a SPIKE
        % via a poisson process with timestep DT
            noisyRate = rate * (1 + 0.1*randn());
            spike = (rand() < (dt * noisyRate));
        end
        
    end
end