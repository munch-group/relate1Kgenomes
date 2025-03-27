#!/usr/bin/env python
import json
import sys
import subprocess
import pprint
import re


def find_closing(text, idx):
    starting = text[idx]
    closing = {'{': '}', '[': ']', '(': ')'}[starting]
    stack = []
    for i, c in enumerate(text[idx:]):
        if c == starting:
            stack.append(i)
        if c == closing:
            stack.pop()
            if not stack:
                return i
    
# print(find_closing('[ [1, 1, 1], 2     ]  ', 0))

target = sys.argv[1]
json_info = subprocess.check_output(['gwf', 'info', target]).decode()
info = json.loads(json_info)
s = pprint.pformat(info)
indent = 0
for line in s.split('\n'):
    m = re.search(r"^(\s*)", line)
    if m:
        indent = max(indent, len(m.group(1)))

i = s.index("'spec':") + 7
j = s.index("}}")
# print(s[:i], '"""', end='') 
# spec = s[i:j]
# spec = '\n'.join(re.sub(r"^\s*'\s*", ' '*(indent-4), x[:-1].replace(r"\n", '')) for x in spec.split('\n'))
# print(spec, end='')
# print('"""', s[j:])
# print(indent)

print('\nDEPENDENTS:')

pprint.pprint(info[target]['dependents'])
print('\nOPTIONS:')

pprint.pprint(info[target]['options'])
print('\nINPUTS:')

pprint.pprint(info[target]['inputs'])
print('\nOUTPUTS:')

pprint.pprint(info[target]['outputs'])
print('\nSPEC:\n"""')
#print('===========')
spec = pprint.pformat(info[target]['spec'])[1:-1]#s[i:j]
lines = []
for line in spec.split('\n'):
    line = line.strip()
    line = line[1:-1]
    line = line.strip()
    line = line.replace('\\n', '\n')
    line = line.lstrip()
    if not line:
        continue
    line = ' ' + line
    lines.append(line)
print(''.join(lines))
print('"""')

#     line = line.strip()
#     # line = line[:-1].replace(r"\n", '')
#     # line.replace(r"\\n", '')
# #    line = re.sub(r"^\s*'\s*", ' '*, line)
#     print(line, end=endl)


# target = sys.argv[1]
# json_info = subprocess.check_output(['gwf', 'info', target]).decode()
# info = json.loads(json_info)
# print(repr(info))
# for line in repr(info).split('\n'):
#     print('xx')
#     if 'spec' in line:
#         break
#     # print(line)


# target = sys.argv[1]
# json_info = subprocess.check_output(['gwf', 'info', target]).decode()
# info = json.loads(json_info)
# s = pprint.pformat(info)
# s = pprint.pformat(info)
# i = s.index("'spec':")
# print(s[:i], end='')
# for line in s[i:].split('\n'):
#     m1 = re.search(r"^(\s*)'(.*)'\s*\Z", line)
#     m2 = re.search(r"^(\s*)'(.*)'}}\Z", line)
#     if m1:
#         line = m1.group(1) + m1.group(2)
#     elif m2:
#         line = f'{m2.group(1)} {m2.group(2)}"""'
#     if line.endswith('\\n'):
#         line = line[:-2]
#     else:
#         line = line + '\\'
# #    line = line.replace('\\n', '')


#     print(line)
#     # if 'spec' in line:
#     #     print(line)
#     #     break
# # re.sub("'spec': ['", '"""', s)
# # print(re.search(r'(?<=spec: )\n', s).groups())



# spec = pprint.pformat(info[target]['spec'].split('\n'))

# print(spec)

#info[target]['spec'] = info[target]['spec'].split('\n')
#pprint.pprint(info)
# print(repr(info[target]['spec'].replace('\\\\n', '\\n')))  #repr(info[target]['spec'])

#print(info)
# with open(f'.gwf/logs/{sys.argv[1]}.stdout') as f:
#     print(json.loads(f.read()))

#    print(json.load(f))
