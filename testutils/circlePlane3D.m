function [X,Y,Z] = circlePlane3D(C, N, R, n_, AL)


N = N ./ repmat(sqrt(dot(N,N,2)),1,3);

U = [zeros(n_,1) ones(n_,1) -N(:,2)./N(:,3)];
U = U ./ repmat(sqrt(dot(U,U,2)),1,3);

V = cross(N,U);

AL = AL(:)';
AL_ = repmat(AL,n_,1);

X = cell(3,1);
R_ = repmat(R,1,length(AL));
for i = 1: 3
    C_ = repmat(C(:,i),1,length(AL));
    U_ = repmat(U(:,i),1,length(AL));
    V_ = repmat(V(:,i),1,length(AL));
    X{i} = C_ + R_.*cos(AL_).*U_ + R_.*sin(AL_).*V_;
end

Y = X{2};
Z = X{3};
X = X{1};


