from pymol import cmd, stored

filename = 'C:\\Python27\\Lib\\site-packages\\pmg_tk\\startup\\MitchellSheet3.tsv'
cmd.delete('all')
#cmd.window('hide')
stored.mystuff = {'aspot':[]}
print 'hello'
pdbsfile = open(filename,'r')
resultsfile = open('results.tsv','w+')
spltln = pdbsfile.readline().split('\t')
while len(spltln)>1:
    print 'hello'
    resultsfile.write('\t'.join(spltln))
    pdb1name = spltln[0]
    resisAndPosis = spltln[2:]
    resis = []
    posis = []
    for e in resisAndPosis:
        if e=='': break
        spliteded = e.split(',')
        resis.append(spliteded[0])
        posis.append(spliteded[1])
    pdb1 = cmd.fetch(pdb1name,async=0)
    resultsfile.write(pdbsfile.readline())
    spltln = pdbsfile.readline().split('\t')
    while spltln[2]!='':
        pdb2name = spltln[2].split()[0]
        if pdb1name==pdb2name or len(pdb2name)!=4:
            spltln = pdbsfile.readline().split('\t')
            continue
        spltwrtln = spltln[:3]
        print pdb2name
        pdb2 = cmd.fetch(pdb2name,async=0)
        if not pdb2name in cmd.get_names():
            spltln = pdbsfile.readline().split('\t')
            continue
        alignment = cmd.align(pdb1name,pdb2name)
        spltwrtln.append("%.3f"%(alignment[0]))
        slctnlst = []
        matches = []
        for i in range(len(resis)):
            disresi = 'a' + str(i)
            slctnlst.append(disresi)
            #cmd.select(disresi,'br. resi %s in %s a. 1 & br. resn %s in %s'%(str(posis[i]),pdb1name,str(resis[i]),pdb2name))
            cmd.select(disresi,'resi %s in %s'%(posis[i],pdb1name))
            cmd.select(disresi,'%s a. 2'%(disresi))
            cmd.select(disresi,'br. %s'%(disresi))
            cmd.select(disresi,'%s & resn %s'%(disresi,resis[i]))
            cmd.select(disresi,'%s in %s'%(disresi,pdb2name))
            cmd.iterate(disresi + " & n. CA","aspot.append(resi)",space=stored.mystuff)
            if cmd.count_atoms(disresi)>0:
                print stored.mystuff['aspot']
                matches.append(i)
                spltwrtln.append('%s,%s'%(resis[i],str(stored.mystuff['aspot'][-1])))
                stored.mystuff['aspot'] = []
            else:
                spltwrtln.append('----------------')
        cmd.select('a',' | '.join(slctnlst))
        blist = []
        for j in matches:
            datresi = 'b' + str(j)
            cmd.select(datresi,'%s a. 2 in %s'%(slctnlst[j],pdb1name))
            blist.append(datresi)
        if blist==[]:
            spltwrtln.insert(4,'-------')
        else:
            cmd.select('c',' | '.join(blist))
            spltwrtln.insert(4,"%.3f"%(cmd.align('c','a')[0]))
        cmd.delete(pdb2name)
        #for f in slctnlst:
            #cmd.delete(f)
        #cmd.delete('a')
        #stored.mystuff['aspot'] = []
        resultsfile.write('\t'.join(spltwrtln) + '\n')
        resultsfile.flush()
        spltln = pdbsfile.readline().split('\t')
    cmd.delete(pdb1name)
    #cmd.delete('b')
    spltln = pdbsfile.readline().split('\t')
    resultsfile.write('\t\t\t\t\t\t\t\t\n')
#cmd.iterate(pdb1name+" and c. A and n. CA","pdb1.append(resn)",space=stored.mystuff)
#cmd.iterate(pdb2name+" and c. A and n. CA","pdb2.append(resn)",space=stored.mystuff)
print '%s %s'%(pdb1name,pdb2name)
pdbsfile.close()
resultsfile.close()
cmd.window('show')
    
#alignFrom('MitchellSheet.tsv')
