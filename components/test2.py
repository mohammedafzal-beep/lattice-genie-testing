dict_key_map = {}
for i in range(1,64):
    if i < 15:
        dict_key_map[i] = [0, i-1]
    elif i < 29:
        dict_key_map[i] = [1, i-15]
    elif i < 36:
        dict_key_map[i] = [2, i-29]
    elif i < 50:
        dict_key_map[i] = [3, i-36]
    elif i < 64:
        dict_key_map[i] = [4, i-50]
print(dict_key_map)
    