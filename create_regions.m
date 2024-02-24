

function [table, cas] = create_regions(A, B, V, N, show)

	if (nargin < 5) || isempty(show)
		show = 0;
	end

	maxIterations = 1000;
	
	Y = (ones(length(A),length(B)).*B);
	X = (ones(length(B),length(A)).*A)';
	
	[finalCentroids, finalAssignments] = kMeans([X(:),Y(:),V(:)], N, maxIterations);
	
	%%
	
	%cdata = reshape(finalAssignments(:), [length(V),length(T)])
	
	finas = finalAssignments(:);
	idx = 1;
	i=1;
	addext = 0;
	cas = finalCentroids;
	
	if finas(i)>1
		store = finas(i);
		addext = addext+1;
		cas = [cas; cas(idx,:)];
		for j=1:length(finas)
			if finas(j)==1
				finas(j)=N+addext;
			end
		end
		for j=i:length(finas)
			if finas(j)==store
				finas(j)=1;
			end
		end
		cas(idx,:) = cas(store,:);
	end
	
	for i=2:length(finas)
		if finas(i)>idx
			if finas(i)==idx+1
				idx = finas(i);
			else
				store = finas(i);
				addext = addext+1;
				idx = idx+1;
				cas = [cas; cas(idx,:)];
				for j=1:length(finas)
					if finas(j)==idx
						finas(j)=N+addext;
					end
				end
				for j=i:length(finas)
					if finas(j)==store
						finas(j)=idx;
					end
				end	
				cas(idx,:) = cas(store,:);
			end
		end
	end
	
	table = reshape(finas(:), [length(A),length(B)]);
	cas = cas(1:N,:);

	if show
		% Assign color based on the assignments
		colormap(parula(max(finalAssignments(:))));
	
		% Set the CData property to color the surface based on assignments
		set(gca, 'CLim', [min(finas(:)), max(finas(:))]); % Adjust color limits
		set(surf(X,Y,V), 'CData', table, 'FaceColor', 'interp');
		
		% Add colorbar for reference
		colorbar;
		
		% Display final centroids and assignments
		disp('Final Centroids:');
		disp(finalCentroids);
	
	%disp('Final Assignments:');
	%disp(finalAssignments);
	end
end

%%

function [centroids, assignments] = kMeans(points, k, maxIterations)
    numPoints = size(points, 1);
    numDimensions = size(points, 2);
    
    % Initialize centroids randomly from the data points
    centroids = points(randperm(numPoints, k), :);
    
    for iter = 1:maxIterations
        % Assign each point to the nearest centroid
        distances = pdist2(points(:,numDimensions), centroids(:,numDimensions));
        [~, assignments] = min(distances, [], 2);
        
        % Update centroids
        for i = 1:k
            clusterPoints = points(assignments == i, :);
            if ~isempty(clusterPoints)
                centroids(i, :) = mean(clusterPoints);
            end
        end
    end
end