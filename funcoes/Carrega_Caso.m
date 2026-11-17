function mpc = Carrega_Caso(cfg)
%CARREGA_CASO  Carrega o sistema definido em cfg.caso no formato do codigo.
%
%   mpc = Carrega_Caso(cfg)
%
%   cfg.caso: 'case4gs'                -> 4 barras (Grainger e Stevenson)
%             'case118'                -> IEEE 118 barras (MATPOWER)
%             'case118_base'           -> IEEE 118 barras (PGLib-OPF) no
%                                         ponto base de Gera_Base118
%             'case_ACTIVSg2000_hNNNN' -> ACTIVSg2000 na hora NNNN da serie
%                                         de 2016 (Gera_BaseACTIVSg2000)
%
%   As areas sao dadas por cfg.mapaAreas e os GCIs por cfg.combGCI e
%   cfg.GCIextra (ver Converte_MATPOWER e CIA_Dados).



    %% a) Leitura e Conversao do Caso
    switch lower(cfg.caso)
        case 'case4gs'
            mpc = Converte_MATPOWER(case4gs, cfg.mapaAreas, {});

        case 'case118'
            mpc = Converte_MATPOWER(case118, cfg.mapaAreas, {});

        case 'case118_base'
            mpc = Converte_MATPOWER(case118_base, cfg.mapaAreas, {});

        otherwise
            if strncmpi(cfg.caso, 'case_ACTIVSg2000_h', 18)
                mpc = Converte_MATPOWER(feval(cfg.caso), cfg.mapaAreas, cfg.combGCI);
            else
                error('Carrega_Caso: caso "%s" desconhecido.', cfg.caso);
            end
    end
end
