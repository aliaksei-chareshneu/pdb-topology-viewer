// TODO: Convert cartesian to Y-reversed coordinates
// TODO: Fix tsc errors
// TODO: Copy updated function to index.ts and continue development of this function in index.ts
// TODO: Implement rotation
// TODO: Check if residues from some of 5 used APIs correspond to what is in 2DProts
// TODO: Write better description
// Function converts 2DProts JSON to "PDBe-topology-API-like" JSON suitable for drawing SSEs via modified PDB Topology Component
function convert2DProtsJSONtoTopologyAPIJSON (inputJson, entryID, chainID) {
	// TODO: try different for both if something goes wrong
	const MINORAXIS = 4;
	const CONVEXITY = 2;
	
	// Coordinates of upper right and lower left corners of "canvas"
	const upperRight = {
		'x': inputJson.metadata['upper_right'][0],
		'y': inputJson.metadata['upper_right'][1]
	};
	
	const lowerLeft = {
		'x': inputJson.metadata['lower_left'][0],
		'y': inputJson.metadata['lower_left'][1]
	};
	
	// const outputJSON = {
		// entryID: {
			// TODO: check if it should be determined in some way (e.g. if chainID = A, entityID = 1, if B => 2, etc.)
			// '1': {
				// chainID: {
					// 'helices': [],
					// 'coils': [],
					// 'strands': [],
					// 'terms': [],
					// 'extents': [],
				// }
			// }
		// }
	// };
	
	const outputJSON = {};
	// TODO: check if entityId (i.e. '1') should be determined in some way (e.g. if chainID = A, entityID = 1, if B => 2, etc.)
	outputJSON[entryID] = {'1': {}};
	outputJSON[entryID]['1'][chainID] = {
		'helices': [],
		'coils': [],
		'strands': [],
		'terms': [],
		'extents': [],
	};
	
	// const helices = outputJSON.entryID['1'].chainID.helices;
	
	const inputSSEs = Object.entries(inputJson.sses);
	const inputHelices = inputSSEs.filter(sse => sse[0].charAt(0) === 'H');
	// console.log(inputHelices);
	outputJSON[entryID]['1'][chainID].helices = inputHelices.map(helix => {
		console.log(helix);
		const center = {
			'x': helix[1].layout[0],
			'y': helix[1].layout[1],
		};
		
		return {
			'start': Number(helix[1].residues[0]),
			'stop': Number(helix[1].residues[1]),
			'majoraxis': Number(helix[1].size),
			'minoraxis': MINORAXIS,
			// We can try 2 points here, and let the existing code of topology component convert it
			// But better to let the code commented, and supply 6 points from the very beginning
			'path': [
				center.x + (MINORAXIS/2), center.y + (helix[1].size/2) - CONVEXITY,
				center.x,				  center.y + (helix[1].size/2),
				center.x - (MINORAXIS/2), center.y + (helix[1].size/2) - CONVEXITY,
				center.x - (MINORAXIS/2), center.y - (helix[1].size/2) + CONVEXITY,
				center.x,				  center.y - (helix[1].size/2),
				center.x + (MINORAXIS/2), center.y - (helix[1].size/2) + CONVEXITY
			],
			'color': helix[1].color,
			'2dprotsSSEId': helix[0],
		}
	});
	
	return outputJSON;
	
	// for (const [key, value] of inputHelices)
}