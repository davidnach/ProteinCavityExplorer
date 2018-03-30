const express = require('express')
const router = express.Router();


router.post('/createExperiment',(req,res) => {   
    console.log(req.body);
    //pid = req.body[0].pid;
    info = [{pid: '1ake', id : '12'}];
	res.send(JSON.stringify(info));
});





module.exports = router;
