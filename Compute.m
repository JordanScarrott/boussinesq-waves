classdef Compute
    properties
        Property1
    end
    
    methods(Static)
        % Takes the required arrays (n,U,V,h) and constant params 
        % (p = [a1 a2 b1 b2 g, dx, dy]) required to compute E
        
        % Test String
        % First Run patel2020_v[#] then enter following
        % compute_E(n(:,:,1), u(i).u, u(i).v, h, [a1 a2 b1 b2 g dx dy])
        function E_new = E(ni, ui, vi, h, p)
            a1i = p(1);
            a2i = p(2);
            b1i = p(3);
            b2i = p(4);
            gi = p(5);
            dxi = p(6);
            dyi = p(7);
            
            % Setting up sub-components pf E
            uxx = doublederiv2d(ui, 1, dxi);
            vyy = doublederiv2d(vi, 2, dyi);
            uxy = deriv2d(deriv2d(ui, 1, dxi), 2, dyi);
            vxy = deriv2d(deriv2d(vi, 1, dxi), 2, dyi);
                    
            huxx = doublederiv2d(h .* ui, 1, dxi);
            hvyy = doublederiv2d(h .* vi, 2, dyi);
            huxy = deriv2d(deriv2d(h .* ui, 1, dxi), 2, dyi);
            hvxy = deriv2d(deriv2d(h .* vi, 1, dxi), 2, dyi);
        
            
            % Computing components of E
            E_1 = deriv2d((h + ni) .* ui, 1, dxi);
            E_2 = deriv2d((h + ni) .* vi, 2, dyi);
            E_3 = a1i * h.^3 .* (uxx + vxy);
            E_4 = a2i * h.^2 .* (huxx + hvxy);
            E_34 = deriv2d(E_3 + E_4, 1, dxi);
            E_5 = a1i * h.^3 .* (vyy + uxy);
            E_6 = a2i * h.^2 .* (hvyy + huxy);
            E_56 = deriv2d(E_5 + E_6, 2, dyi);
            
            
            % Adding together final components of E and sizing the arrays correctly
            E_new = -E_1 - E_2 - E_34 - E_56;
        end

        function F_new = F(ni, ui, vi, g, dxi, dyi)
            ux = deriv2d(ui, 1, dxi);
            uy = deriv2d(ui, 2, dyi);
            nix = deriv2d(ni, 1, dxi);
        
            F_new = -g * nix - (ui .* ux + vi .* uy);
        end

        function F1_new = F1(h, vi, dxi, dyi, b1, b2)
            vx = deriv2d(vi, 1, dxi);
            vxy = deriv2d(vx, 2, dyi);
            
            hvx = deriv2d(h .* vi, 1, dxi);
            hvxy = deriv2d(hvx, 2, dyi);
            
            F1_new = -h .* (b1*(h .* vxy) + b2*hvxy);
        end

        function G_new = G(ni, ui, vi, g, dxi, dyi)
            vx = deriv2d(vi, 1, dxi);
            vy = deriv2d(vi, 2, dyi);
            niy = deriv2d(ni, 2, dyi);
        
            G_new = -g * niy - (vi .* vy + ui .* vx);
        end

        function G1_new = G1(h, ui, dxi, dyi, b1, b2)
            ux = deriv2d(ui, 1, dxi);
            uxy = deriv2d(ux, 2, dyi);
            
            hux = deriv2d(h .* ui, 1, dxi);
            huxy = deriv2d(hux, 2, dyi);
            
            G1_new = -h .* (b1*(h .* uxy) + b2*huxy);
        end

        function nt_new = nt(ni,ui,h,za_coefficient,dx,dy)
            dims_ui = size(ui);
                
            za = za_coefficient * h;
            
            % Divergence of term 1
            diverg1 = (h+ni).*ui.u + (h+ni).*ui.v;
            
            % Coefficients for terms 2 and 3
            c2 = (za.^2/2 - h.^2/6);
            c3 = (za + h/2);
            
            % Term 2 - h * grad (grad dot u)
            term2 = Vector_Field([dims_ui(1), dims_ui(2)]);
            ux = deriv2d(ui.u,1,dx);
            vy = deriv2d(ui.v,2,dy);
            term2.u = c2 .* h .* deriv2d(ux + vy,1,dx);
            term2.v = c2 .* h .* deriv2d(ux + vy,2,dy);
            
            
            % Term 3 - h * grad (grad dot (hu))
            term3 = Vector_Field([dims_ui(1), dims_ui(2)]);
            hux = deriv2d(h.*ui.u, 1, dx);
            hvy = deriv2d(h.*ui.v, 2, dy);
            term3.u = c3 .* h .* deriv2d(hux, 1, dx);
            term3.v = c3 .* h .* deriv2d(hvy, 2, dy);
            
            % Divergence of (term2 + term3)
            diverg2 = deriv2d(term2.u + term3.u,1,dx) + deriv2d(term2.v + term3.v,2,dy);
            
            % Final Answer
            nt_new = -diverg1 - diverg2;
        end

        function U_new = U(u, h, b1, b2, dx)
            U_new = u + h .* (b1 .* h .* doublederiv2d(u,1,dx) + b2 .* doublederiv2d(h .* u,1,dx));
        end

        function V_new = cV(v, h, b1, b2, dy)
            V_new = v + h .* (b1 .* h .* doublederiv2d(v,2,dy) + b2 .* doublederiv2d(h .* v,2,dy));
        end
    end
end

