import v from 'vector-math'

const { createVectorObj, crossProduct, dotProduct, subVector, unitVector } = v


export default function parse_mmCIF(text) {

    const PDB = text.split('\n')
    const ATOMS = text.split('\n').filter((line) => line.startsWith('ATOM')).map(array => array.split(' ').filter(element => element !== ''))
    const customLabels = ['', 'id', 'atom', 'atom_type', '', 'residue', 'chain', 'entity_index', 'residue_index', '', 'x', 'y', 'z', 'occupancy', 'isotropic_temperature_factor', 'formal_charge', 'author_residue_index', 'author_residue', 'author_chain', 'author_atom_type', '' ]    
    
    let object = {}
    const newAtoms = ATOMS.map((array) => {
        let atom_info = {}
        let author_entries = {}
        array.forEach((entry, index) => {
            const label = customLabels[index]
            if(label !== '') {
                if(label.includes('author')) {
                    if( label.includes('index')) {
                        const int = parseInt(entry)
                        author_entries = {...author_entries, [label.replace('author_', '')]: int}
                    }
                    else {
                        author_entries = {...author_entries, [label.replace('author_', '')]: entry}
                    }
                    atom_info = {...atom_info, author_entries}
                    return
                }
                if(label === 'id' || label.includes('index') || label === 'formal_charge') {
                    const int = parseInt(entry)
                    atom_info = ({...atom_info, [label]: int})
                    return
                }
                if(label === 'isotropic_temperature_factor' || label === 'x' || label === 'y' || label === 'z' || label === 'occupancy') {
                    const float = parseFloat(entry)
                    atom_info = ({...atom_info, [label]: float})
                    return
                }
                else {
                    atom_info = ({...atom_info, [label]: entry})
                }
            }
        })
        return (
            atom_info
        )
    })
    const entityStart = PDB.indexOf('_entity.details ')
    const chainData = []

    for(let i=entityStart; i <= 1000; i++) {
        if(PDB[i] === '# ') {
            break
        } else {
            chainData.push(PDB[i])
        }
    }

    chainData.shift()

    const chainInfo = chainData.map((line) => replaceSpaces(line, 0)).map((line) => replaceQuotes(line, 0)).map((line) => line.split(' ')).map((array) => array.filter((element) => element !== ''))

    function replaceSpaces(line, count) {
        if(count === line.length-1) {
            return (
                line
            )
        }
        const newLine = line.replace(/(?<='\w+)\s/, '_')
        count++
        return replaceSpaces(newLine, count)
    }

    function replaceQuotes(line, count) {
        if(count >= 2) {
            return (
                line
            )
        }
        const newLine = line.replace(/['"`]/, '')
        count++
        return replaceQuotes(newLine, count)
    }


    const chainArray = chainInfo.map((entry) => {
        let chainObj = {}
        entry.forEach((entry, index) => {
            switch (index) {
                case 1:
                    chainObj = {...chainObj, type: entry }
                    break;
                case 3:
                    chainObj = {...chainObj, name: entry }
                    break;
                case 4:
                    chainObj = {...chainObj, molecular_weight: parseFloat(entry) }
                    break;
                case 5:
                    chainObj = {...chainObj, quantity: parseInt(entry) }
                    break;
            
                default:
                    break;
            }
        })
        return chainObj
    })


    object = {...object, chain_info: chainArray}
    
    let chainCount = 0
    object.chain_info.filter((entry) => entry.type === 'polymer').forEach((entry) => chainCount += entry.quantity)

    const chains = []
    for(let i=0; i < chainCount; i++) {
        chains.push([])
    }

    const chainLabels = newAtoms.map((atom) => atom.chain).filter((chainID, index, array) => chainID !== array[index+1])
    const chainAtoms = chains.map((chain, index) => newAtoms.filter((atom) => atom.chain === chainLabels[index]))

    let chainObject = {}
    chainAtoms.forEach((chain, index) => {
        chainObject = {...chainObject, [chainLabels[index]]: chain } })

    object = {...object, chains: chainObject}
    object = {...object, atoms: newAtoms}


    const backbones = chainAtoms.map((chain) => chain.filter((residue) => residue.atom_type === 'CA' || residue.atom_type === 'C' || residue.atom_type === 'N' ))
    
    let backboneObject = {}
    backbones.forEach((chain, index) => {
        backboneObject = {...backboneObject, [chainLabels[index]]: chain }
    })

    object = {...object, backbones: backboneObject }


    // torsion angles

    const residues = backbones.map((chain) => {
        return (
            chain.filter((atom, index, array) => {
                if(index !== array.length-1) {
                return (
                    atom.residue !== array[index+1].residue
                )
            }
        }).map((atom) => atom.residue)
        )
    })

    const residueAtoms = residues.map((chain, index) => chain.map((residue, i) => backbones[index].slice(i*3, i*3+3) ))
    const residueAtomsWLabels = residues.map((chain, index) => chain.map((residue, i) => {const resObj = { [residue]: backbones[index].slice(i*3, i*3+3) }; return resObj}))



    const torsionAngles = residueAtoms.map((chain, i) => {
        return (
            chain.map((residue, index, array) => {
                if(index !== 0 && index !== array.length-1) {
                    const vectors = { 
                        Ni: createVectorObj([residue[0].x, residue[0].y, residue[0].z]),

                        Cix: createVectorObj([array[index-1][2].x, array[index-1][2].y, array[index-1][2].z]), 
                        Cia: createVectorObj([residue[1].x, residue[1].y, residue[1].z]),

                        Nii: createVectorObj([array[index+1][0].x, array[index+1][0].y, array[index+1][0].z]), 
                        Ci: createVectorObj([residue[2].x, residue[2].y, residue[2].z]) 
                    }

                    
                    const phiPlanes = [[vectors.Cix, vectors.Ni, vectors.Cia], [vectors.Ni, vectors.Cia, vectors.Ci]]
                    const psiPlanes = [[vectors.Ni, vectors.Cia, vectors.Ci], [vectors.Cia, vectors.Ci, vectors.Nii]]
                    
                    const phiNormals = phiPlanes.map((plane) => {
                        const U = v.subVector(plane[0], plane[1])
                        const W = v.subVector(plane[0], plane[2])
                        const V = v.crossProduct(U, W)
                        const Vi = v.unitVector(V)

                        return Vi
                    })

                    const psiNormals = psiPlanes.map((plane) => {
                        const U = v.subVector(plane[0], plane[1])
                        const W = v.subVector(plane[0], plane[2])
                        const V = v.crossProduct(U, W)
                        const Vi = v.unitVector(V)
                        
                        return Vi
                    })

                    const phiSign = Math.sign((phiNormals[1].i - phiNormals[0].i)*(phiNormals[1].j - phiNormals[0].j)*(phiNormals[1].k - phiNormals[0].k))
                    const psiSign = Math.sign((psiNormals[1].i - psiNormals[0].i)*(psiNormals[1].j - psiNormals[0].j)*(psiNormals[1].k - psiNormals[0].k))

                    const phi = Math.acos(v.dotProduct(phiNormals[0], phiNormals[1])) * 180/Math.PI*phiSign
                    const psi = Math.acos(v.dotProduct(psiNormals[0], psiNormals[1])) * 180/Math.PI*psiSign

                    const angles = {residue: residue[0].residue, phi: phi, psi: psi}

                    return angles
                    
                } else return null
            })
        )
    })

    let torsionObject = {}
    residues.forEach((chain, index) => { torsionObject = {...torsionObject, [chainLabels[index]]: chain.map((residue, i, array) => { if(i !== 0 && i !== array.length-1) {  return torsionAngles[index][i] } else { const blank =  { residue: residue, phi: null, psi: null}; return blank } })} } )

    object = { ...object, torsion_angles: torsionObject }
    
    return object

}
