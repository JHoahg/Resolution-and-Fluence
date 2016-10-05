%Author JH
% class for propagator, reduces comp over head if used correctly
% example:
% prop = Propagator(1e-3, 1e-3, 2048, 2048); % only once
% % opt 5th parameter for size of padding
% bla = phantom(2048);
% bla2 = prop.propTF(bla);
% imagesc(abs(bla2));

classdef PropagatorGPU
    
    properties (SetAccess = 'private')
        fresnelx;
        fresnely;
        f; %enlargement factor
        nx;
        ny;
        H;
        gpu_kernel;
    end
    
    
    methods
        function obj = PropagatorGPU(fresnelx,fresnely, nx, ny, f, do_whats_necessary)
            obj.fresnelx = fresnelx;
            obj.fresnely = fresnely;
            
            if (nargin < 5)
                obj.f = 0;
            else
                obj.f = f;
            end
            
            obj.nx = nx;
            obj.ny = ny;
            
            
            
            obj = init(obj);
            
        end %end constructor
        
        function obj = init(obj)
            
            obj.gpu_kernel = parallel.gpu.CUDAKernel( 'prop_kernel.ptx', ...
  'prop_kernel.cu' );
            
%             warning('propagator reinit')
            % Sampling criterion for phase chirp.
            px = abs(1/(obj.nx*obj.fresnelx));
            py = abs(1/(obj.ny*obj.fresnely));
            if (obj.f == 0)
                obj.f = max(px,py);
                if (obj.f < 1)
                    obj.f = 1;
                end
                %     do_whats_necessary = 0;
            end
            
            
            if (max(obj.nx,obj.ny)*obj.f > 100000)
                yn = input(sprintf('Enlarged image would be very large (%ix%i Pixels).\nReally use this factor f=%f? (N) ',round(obj.nx*obj.f),round(obj.ny*obj.f),obj.f),'s');
                if ~(strcmpi(yn,'yes') || strcmpi(yn,'y'))
                    return;
                end
            end
            
            if obj.f > 1
                fprintf('Probe-FOV is enlarged by factor %3.2f for propagation.\n',obj.f);
                
                % enlarge array by these sizes
                ny_l = round((ceil(obj.f)*obj.ny));
                nx_l = round((ceil(obj.f)*obj.nx));
%                 im_l = padarray(im,[pady padx],'replicate');
                %TODO: we should test here for odd numbers
                
                % new coordinate system
%                 [ny_l, nx_l] = size(im_l);
%                 dqy = 2*pi/(ny_l);
%                 dqx = 2*pi/(nx_l);
                
                
%                 [Qx,Qy] = meshgrid(((1:nx_l)-floor(nx_l/2)-1)*dqx,((1:ny_l)-floor(ny_l/2)-1)*dqy);
%                 % paraxial propagator
%                 kappa_x = -1/2*ifftshift(Qx.^2);
%                 kappa_y = -1/2*ifftshift(Qy.^2);
                
                % propagation kernel
%                 obj.H = exp(1i*kappa_x/(2*pi*obj.fresnelx)+1i*kappa_y/(2*pi*obj.fresnely));
%                 obj.H=fftshift(obj.H); %shift trans func
                obj.nx = nx_l;
                obj.ny = ny_l;
                
            else
                
%                 dqx = 2*pi/(obj.nx);
%                 dqy = 2*pi/(obj.ny);
%                 
%                 [Qx,Qy] = meshgrid(((1:obj.nx)-floor(obj.nx/2)-1)*dqx,((1:obj.ny)-floor(obj.ny/2)-1)*dqy);
%                 kappa_x = -1/2*ifftshift(Qx.^2);
%                 kappa_y = -1/2*ifftshift(Qy.^2);
                
                % propagation kernel
%                 obj.H = exp(1i*kappa_x/(2*pi*obj.fresnelx)+1i*kappa_y/(2*pi*obj.fresnely));
%                 obj.H=fftshift(obj.H); %shift trans func
                
            end 
            
%             obj.H = gpuArray(obj.H);
        end
        
        
        function[u1]=propTF(obj, u1)  
          if obj.fresnelx < Inf
            if obj.f > 1
                %fprintf('Probe-FOV is enlarged by factor %3.2f for propagation.\n',obj.f);
                
                [M,N]=size(u1); %get input field array size
                
                % enlarge array by these sizes
%                 pady = round((obj.ny*f - M)/2);
%                 padx = round((obj.nx*f - N)/2);
                %TODO: modify padsmoothly to use pad to size...
                pady = round((ceil(obj.f)*M -M)/2);
                padx = round((ceil(obj.f)*N -N)/2);
                u1 = padarray(u1,[pady padx],'replicate');
                
%                 if (M ~= obj.ny || N ~= obj.nx )
%                     obj.nx = N;
%                     obj.ny = M;
%                     obj = init(obj);
%                     warning('propagator reinit')
%                 end
            end
            
%             u1 = fftshift(ifft2(fft2(ifftshift(u1)).*obj.H));
%             U2=fft2(fftshift(u1)); %shift, fft src field
%             U2=obj.H .* U2; %multiply
%             U2=ifftshift(ifft2(U2)); %inv fft, center obs field
            
%     configure execution grid
            block_size = 512;
            obj.gpu_kernel.ThreadBlockSize = [block_size, 1, 1];
            obj.gpu_kernel.GridSize = [ceil(obj.ny * obj.nx  / block_size)];
            
            u1 = fft2(u1);
            u1 = feval(obj.gpu_kernel, u1, u1, obj.fresnelx, obj.ny, obj.nx,...
                       obj.ny*obj.nx);
            u1 = ifft2(u1);
            
            if obj.f>1
                % decrease field of view
                u1 = u1(pady+1:end-pady,padx+1:end-padx);
            end
          end
        end %proptf end
        
        function[in]=kernel(obj, in)
           in = obj.H; 
        end
    end %methods end
    
end