__author__ = 'curtisj'

import shlex

DEBUG = False

ignore_list = ['>', '<', '!=', '==', '>=',
               '<=', '=', 'and', 'not', 'or', '(', ')']
word_list = ['and', 'not', 'or']


def parse_basis(basis):

    kd = {}
    kd['atom'] = ['atom[i] == ', 'atom[i] ']
    kd['index'] = ['index[i] == ', 'index[i] ']
    kd['name'] = ['name[i] == ', 'name[i] ']
    kd['resname'] = ['resname[i] == ', 'resname[i] ']
    kd['chain'] = ['chain[i] == ', 'chain[i] ']
    kd['resid'] = ['resid[i] == ', 'resid[i] ']
    kd['occupancy'] = ['occupancy[i] == ', 'occupancy[i] ']
    kd['beta'] = ['beta[i] == ', 'beta[i] ']
    kd['segname'] = ['segname[i] == ', 'segname[i] ']
    kd['element'] = ['element[i] == ', 'element[i] ']
    kd['charge'] = ['charge[i] == ', 'charge[i] ']
    kd['moltype'] = ['moltype[i] == ', 'moltype[i] ']

    if DEBUG:
        print 'original_basis = ', basis, '\n'

    lexer = shlex.shlex(basis)
    original_tokenlist = []
    for token in lexer:
        original_tokenlist.append(str(token))
    if DEBUG:
        print 'original_tokenlist = ', original_tokenlist

    tokenlist = []
    number_of_tokens = len(original_tokenlist)
    if DEBUG:
        print 'len(original_tokenlist) = ', len(original_tokenlist)

    i = 0
    while (i < number_of_tokens):
        this_token = original_tokenlist[i]
        if DEBUG:
            print 'this_token = ', this_token,
        if i < number_of_tokens - 1:
            next_token = original_tokenlist[i + 1]
            if DEBUG:
                print 'next_token = ', next_token
            if(this_token in ignore_list):
                if(next_token in ignore_list):
                    if(this_token in word_list or next_token in word_list):
                        tokenlist.append(this_token + ' ' + next_token)
                    else:
                        tokenlist.append(this_token + next_token)
                    i += 1
                elif(this_token == '='):
                    tokenlist.append(this_token + '=')
                    i += 1
                else:
                    tokenlist.append(this_token)

            elif(this_token not in kd and not this_token.isdigit()):
                tokenlist.append('"' + this_token + '"')
            else:
                tokenlist.append(this_token)
        else:
            if(this_token not in kd and not this_token.isdigit()):

                if(this_token[0] != '"' and this_token not in ignore_list):
                    tokenlist.append('"' + this_token + '"')
                else:
                    tokenlist.append(this_token)
            else:
                tokenlist.append(this_token)
        i += 1

    if DEBUG:
        print 'tokenlist = ', tokenlist

    number_of_tokens = len(tokenlist)
    new_basis = ''

    for i in xrange(number_of_tokens):
        this_word = tokenlist[i]
        try:
            next_word = tokenlist[i + 1]
            nwe = True
        except:
            nwe = False

        if this_word in kd:
            if DEBUG:
                print 'tw,nw = ', this_word, next_word
            if nwe:
                if next_word not in ignore_list:
                    new_basis += kd[this_word][0]
                else:
                    new_basis += kd[this_word][1]

        else:
            new_basis += ' ' + this_word + ' '

    if DEBUG:
        print 'new_basis = ', new_basis, '\n'

    return new_basis

if __name__ == '__main__':
    basis = []
    basis.append('(name CA and name NH) or resid > 43')
    basis.append('(name CA and name NH) or resid > 43 and resid < 57')
    basis.append('segname HC1 and (resid >= 210 and resid <=214)')
    basis.append('segname HC1 and resid < 210')
    basis.append('(resid > 23 and resid < 68) and name "CA"')

    for i in xrange(5):
        print '#####'
        new_basis = parse_basis(basis[i])
        print
        print

    import sassie.sasmol.sasmol as sasmol
    m = sasmol.SasMol(0)
    m.read_pdb('min3.pdb')

    basis = '(resid > 23 and resid < 68) and name "CA"'
    python_basis = parse_basis(basis)

    sub_mol = sasmol.SasMol(0)

    frame = 0
    error, mask = m.get_subset_mask(python_basis)

    if(len(error) > 0):
        print 'error = ', error

    import numpy
    print numpy.sum(mask)

    if(numpy.sum(mask) > 0):
        error = m.copy_molecule_using_mask(sub_mol, mask, frame)

    if(len(error) > 0):
        print 'error = ', error
    else:
        sub_mol.write_pdb('new.pdb', frame, 'w')
