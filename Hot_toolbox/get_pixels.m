function pixels = get_pixels(image)

[row,col] = find(image);
        
pixels = [row,col];