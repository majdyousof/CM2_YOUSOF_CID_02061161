function [bwithout, bwith] = initialiseBridgeModels(E, A, P, L, fy, I)
    b = BRIDGE(E, A, P, L, fy, I);
    bwithout = b.Preprocessor().Assembler();
    bwith = b.Preprocessor_Extra().Assembler();
end
