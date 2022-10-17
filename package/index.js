import v from 'vector-math'

const { createVectorObj, crossProduct, dotProduct, subVector, unitVector } = v

export default function parse_mmCIF () {

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


    // seq with holes begin
    const sequenceStart = PDB.indexOf('_entity_poly_seq.hetero ')
    const seq = []

    for(let i=sequenceStart; i <= 100000; i++) {
        if(PDB[i] === '# ') {
            break
        } else {
            seq.push(PDB[i])
        }
    }


    const polymers = chainInfo.filter((chain) => chain[1] === 'polymer')
    const sequences = polymers.map((chain, index) => seq.filter((line) => line.startsWith(index+1)).map((line) => line.split(' ').filter((entry) => entry !== '')[2]))

    const seqChains = chainLabels.map((chain, index) => sequences[chainAtoms[index][0].entity_index-1])

    const seqResidues = chainAtoms.map((chain, index) => {
            const residueArray = [[]]
            let n = 0
            chain.forEach((atom, index, array) => {
                if(index !== array.length-1) {
                    if(atom.residue_index === array[index+1].residue_index) {
                        residueArray[n].push(atom)
                    } else {
                        residueArray.push([])
                        n++
                    }
                } else { residueArray[n].push(atom) }
            })
        return (
            residueArray
        )
    })


    const chainSequences = seqChains.map((seq, index) => {
        let n = 0
        return (
            seq.map((residue, i) => {            
            if(n <= seqResidues[index].length-1) { 
                    if(residue === seqResidues[index][n][0].residue && i === seqResidues[index][n][0].residue_index-1) {
                        n++
                        return seqResidues[index][n-1]
                    } else {
                        return null
                    }
                } else return null
            })
        )
    })


    const residueBackbones = chainSequences.map((chain, index) => chain.map((residue, i) => {
        if(!residue) {
            return null
        } else return residue.filter((atom) => atom.atom_type === 'N' || atom.atom_type === 'C' || atom.atom_type === 'CA')
    }))

    let chainObject = {}
    chainSequences.forEach((chain, index) => {
        chainObject = {...chainObject, [chainLabels[index]]: chain } })

    let backbonesObject = {}
    residueBackbones.forEach((chain, index) => {
        chainObject = {...chainObject, [chainLabels[index]]: chain } })


    object = {...object, atoms: newAtoms}
    object = {...object, chains: chainObject}
    object = {...object, backbones: backbonesObject }



    // torsion angles

    const torsionAngles = residueBackbones.map((chain, i) => {
        return (
            chain.map((residue, index, array) => {
                    if(residue && array[index-1] && array[index+1]) {
                        const vectors = { 
                            Ni: createVectorObj([residue[0].x, residue[0].y, residue[0].z]),

                            Cix: createVectorObj([array[index-1][2].x, array[index-1][2].y, array[index-1][2].z]), 
                            Cia: createVectorObj([residue[1].x, residue[1].y, residue[1].z]),

                            Nii: createVectorObj([array[index+1][0].x, array[index+1][0].y, array[index+1][0].z]), 
                            Ci: createVectorObj([residue[2].x, residue[2].y, residue[2].z]) 
                    }


                    const phiPlanes = [[vectors.Ni, vectors.Cia, vectors.Ci], [vectors.Cix, vectors.Ni, vectors.Cia]]
                    const psiPlanes = [[vectors.Ni, vectors.Cia, vectors.Ci], [vectors.Cia, vectors.Ci, vectors.Nii]]

                    const phiDirection = v.unitVector(v.subVector(vectors.Cia, vectors.Ni))
                    const psiDirection = v.unitVector(v.subVector(vectors.Cia, vectors.Ci))

                    const phiAlphaCarbonyl = v.unitVector(v.subVector(vectors.Cia, vectors.Ci))
                    const psiAlphaNitrogen = v.unitVector(v.subVector(vectors.Cia, vectors.Ni))

                    
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

                    const phi = transformVectors(phiNormals, phiDirection, phiAlphaCarbonyl)
                    const psi = transformVectors(psiNormals, psiDirection, psiAlphaNitrogen)

                    const angles = {phi: phi, psi: psi}


                    if(i === 0 && (index === 2 || index == 10 || index == 11 || index == 59 || index == 262 )) {
                        transformVectors(psiNormals, psiDirection, psiAlphaNitrogen)
                    }

                    return angles

                } else return null
            })
        )
    })



    let torsionObject = {}
    chainSequences.forEach((chain, index) => { torsionObject = {...torsionObject, [chainLabels[index]]: chain.map((residue, i, array) => { if(residue && array[i-1] && array[i+1]) {  return torsionAngles[index][i] } else { const blank =  { residue: residue, phi: null, psi: null}; return blank } })} } )

    object = { ...object, torsion_angles: torsionObject }



    function multiplyMatrix(inputMatrix, transformMatrix, outputMatrix) {

        outputMatrix[0] = inputMatrix[0] * transformMatrix[0][0] + inputMatrix[1] * transformMatrix[1][0] + inputMatrix[2] * transformMatrix[2][0]
        outputMatrix[1] = inputMatrix[0] * transformMatrix[0][1] + inputMatrix[1] * transformMatrix[1][1] + inputMatrix[2] * transformMatrix[2][1]
        outputMatrix[2] = inputMatrix[0] * transformMatrix[0][2] + inputMatrix[1] * transformMatrix[1][2] + inputMatrix[2] * transformMatrix[2][2]

    }


    function transformVectors(vectors, direction, group) {
                const vectorA = [vectors[0].i, vectors[0].j, vectors[0].k]
                const vectorB = [vectors[1].i, vectors[1].j, vectors[1].k]
                const vectorP = [direction.i, direction.j, direction.k]
                const vectorQ = [group.i, group.j, group.k]

                const thetaZ = Math.atan(-vectorP[0]/vectorP[1])
                const rotationMatrix_Z = [
                        [Math.cos(thetaZ), -Math.sin(thetaZ), 0],
                        [Math.sin(thetaZ), Math.cos(thetaZ), 0],
                        [0, 0, 1]
                    ]

                const intermediateP1 = []
                const intermediateA1 = []
                const intermediateB1 = []
                const intermediateQ1 = []
                multiplyMatrix(vectorP, rotationMatrix_Z, intermediateP1)
                multiplyMatrix(vectorA, rotationMatrix_Z, intermediateA1)
                multiplyMatrix(vectorB, rotationMatrix_Z, intermediateB1)
                multiplyMatrix(vectorQ, rotationMatrix_Z, intermediateQ1)


                const thetaX = Math.atan(intermediateP1[2]/intermediateP1[1])
                const rotationMatrix_X = [
                        [1, 0, 0],
                        [0, Math.cos(thetaX), -Math.sin(thetaX)],
                        [0, Math.sin(thetaX), Math.cos(thetaX)]
                    ]
                        
                const intermediateP2 = []
                const intermediateA2 = []
                const intermediateB2 = []
                const intermediateQ2 = []
                multiplyMatrix(intermediateP1, rotationMatrix_X, intermediateP2)
                multiplyMatrix(intermediateA1, rotationMatrix_X, intermediateA2)
                multiplyMatrix(intermediateB1, rotationMatrix_X, intermediateB2)
                multiplyMatrix(intermediateQ1, rotationMatrix_X, intermediateQ2)


                const thetaY = Math.atan(intermediateQ2[0]/intermediateQ2[2])
                const rotationMatrix_Y = [
                        [Math.cos(thetaY), 0, Math.sin(thetaY)],
                        [0, 1, 0],
                        [-Math.sin(thetaY), 0, Math.cos(thetaY)]
                    ]
                
                const intermediateA3 = []
                const intermediateB3 = []
                const intermediateQ3 = []
                const intermediateP3 = []
                multiplyMatrix(intermediateP2, rotationMatrix_Y, intermediateP3)
                multiplyMatrix(intermediateA2, rotationMatrix_Y, intermediateA3)
                multiplyMatrix(intermediateB2, rotationMatrix_Y, intermediateB3)
                multiplyMatrix(intermediateQ2, rotationMatrix_Y, intermediateQ3)


                const theta180Y = Math.PI
                const rotationMatrix_180Y = [
                    [Math.cos(theta180Y), 0, Math.sin(theta180Y)],
                    [0, 1, 0],
                    [-Math.sin(theta180Y), 0, Math.cos(theta180Y)]
                ]
                const theta180Z = Math.PI
                const rotationMatrix_180Z = [
                    [Math.cos(theta180Z), -Math.sin(theta180Z), 0],
                    [Math.sin(theta180Z), Math.cos(theta180Z), 0],
                    [0, 0, 1]
                ]

                let intermediateA4 = []
                let intermediateB4 = []
                let intermediateQ4 = []
                let intermediateP4 = []
                if(intermediateQ3[2] < 0) {
                    multiplyMatrix(intermediateP3, rotationMatrix_180Y, intermediateP4)
                    multiplyMatrix(intermediateA3, rotationMatrix_180Y, intermediateA4)
                    multiplyMatrix(intermediateB3, rotationMatrix_180Y, intermediateB4)
                    multiplyMatrix(intermediateQ3, rotationMatrix_180Y, intermediateQ4)
                } else {
                    intermediateP4 = intermediateP3
                    intermediateA4 = intermediateA3
                    intermediateB4 = intermediateB3
                    intermediateQ4 = intermediateQ3
                }


                let intermediateA5 = []
                let intermediateB5 = []
                let intermediateQ5 = []
                let intermediateP5 = []
                if(intermediateQ3[1] < 0) {
                    multiplyMatrix(intermediateP4, rotationMatrix_180Y, intermediateP5)
                    multiplyMatrix(intermediateA4, rotationMatrix_180Z, intermediateA5)
                    multiplyMatrix(intermediateB4, rotationMatrix_180Z, intermediateB5)
                    multiplyMatrix(intermediateQ4, rotationMatrix_180Z, intermediateQ5)
                } else {
                    intermediateP5 = intermediateP4
                    intermediateA5 = intermediateA4
                    intermediateB5 = intermediateB4
                    intermediateQ5 = intermediateQ4
                }

                const finalMatrixP = intermediateP5
                const finalMatrixA = intermediateA5
                const finalMatrixB = intermediateB5
                const finalMatrixQ = intermediateQ5



                const localVectors = [v.createVectorObj(finalMatrixA), v.createVectorObj(finalMatrixB)]
                const angle = Math.acos(v.dotProduct(localVectors[0], localVectors[1])) * 180/Math.PI
                const cross = v.crossProduct(localVectors[0], localVectors[1])

                const fixedMatrixA = finalMatrixA.map((float) => parseFloat(float.toFixed(3)))
                const fixedMatrixB = finalMatrixB.map((float) => parseFloat(float.toFixed(3)))
                const fixedVectors = [v.createVectorObj(fixedMatrixA), v.createVectorObj(fixedMatrixB)]
                const fixedAngle = parseFloat(angle.toFixed(3))
                const fixedCross = v.createVectorObj([parseFloat(cross.i.toFixed(3)), parseFloat(cross.j.toFixed(3)), parseFloat(cross.k.toFixed(3))])
                const sign = Math.sign(fixedCross.j)


                return fixedAngle*sign
    }

    return object
}