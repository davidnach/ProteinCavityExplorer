const express = require('express');
const router = express.Router();
var spawn = require("child_process").spawn;





router.post('/createExperiment',  (req,res) => {   
    console.log(req.body);
    var pid = req.body.pid;
    var expID = "";
    var pythonProcess = spawn('python',["pidDirectory.py", pid]);
    pythonProcess.stdout.on('data', function (data){
        // do something with experiment ID (data)
         expID += data;
	 console.log(expID);
         var info = {expId : expID };
         console.log(info);
         res.send(JSON.stringify(info));
    });
});



router.post('/retrieveExperiment', (req,res) => {
 	console.log(req.body);
	var info = {expId : 'experimentid recieved by server'};
	res.send(JSON.stringify(info));
 });	






module.exports = router;
